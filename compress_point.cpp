#include <gmpxx.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Point
{
    mpq_class x;
    mpq_class y;
    bool operator==(const Point& rhs) const {
        return x == rhs.x && y == rhs.y;
    }
};

struct CompressedPoint
{
    mpz_class d; // squarefree signed part of numerator
    mpz_class r; // sqrt(abs(numerator / d))
    mpz_class s; // sqrt(denominator)
};

struct Curve
{
    mpz_class n;
    mpz_class A;
    mpz_class B;

    explicit Curve(const mpz_class &n_)
        : n(n_),
          A(4 * n * n + 12 * n - 3),
          B(32 * (n + 3)) {
    }
};

// ============================================================
// Basic helpers
// ============================================================

static mpz_class absz(const mpz_class &x) {
    return x < 0 ? -x : x;
}

static mpz_class gcd(
    const mpz_class &a,
    const mpz_class &b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}

static bool is_prime(const mpz_class &n) {
    return mpz_probab_prime_p(n.get_mpz_t(), 25) != 0;
}

static mpz_class mod_mul(
    const mpz_class &a,
    const mpz_class &b,
    const mpz_class &m) {
    mpz_class r;
    mpz_mul(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    mpz_mod(r.get_mpz_t(), r.get_mpz_t(), m.get_mpz_t());
    return r;
}

static mpz_class mod_pow(
    const mpz_class &a,
    const mpz_class &e,
    const mpz_class &m) {
    mpz_class r;
    mpz_powm(r.get_mpz_t(), a.get_mpz_t(), e.get_mpz_t(), m.get_mpz_t());
    return r;
}

// ============================================================
// Pollard-rho factorization
// ============================================================

static std::mt19937_64 rng(std::random_device{}());

static mpz_class random_below(
    const mpz_class &n) {
    if (n <= 1)
        return 0;
    mpz_class r;
    const unsigned bits = mpz_sizeinbase(n.get_mpz_t(), 2);
    while (true) {
        mpz_class x = 0;
        unsigned generated = 0;
        while (generated < bits) {
            unsigned long v = rng();
            unsigned take = std::min<unsigned>(64, bits - generated);
            if (take < 64) v &= ((1ul << take) - 1);
            mpz_class part = mpz_class(v);
            mpz_mul_2exp(part.get_mpz_t(), part.get_mpz_t(), generated);
            x += part;
            generated += take;
        }
        x %= n;
        if (x < n) return x;
    }
}

static mpz_class pollard_rho(const mpz_class &n) {
    if (n <= 1) return 1;
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    for (unsigned c = 1;; c++) {
        mpz_class x = random_below(n - 2) + 2;
        mpz_class y = x;
        mpz_class d = 1;
        auto f = [&](const mpz_class &z) -> mpz_class {
            mpz_class result = (mod_mul(z, z, n) + c) % n;
            return result;
        };
        while (d == 1) {
            x = f(x);
            y = f(f(y));
            mpz_class diff = absz(x - y);
            d = gcd(diff, n);
        }
        if (d != n) return d;
    }
}

static void factor_recursive(
    const mpz_class &n,
    std::vector<mpz_class> &factors) {
    if (n == 1) return;
    if (is_prime(n)) {
        factors.push_back(n);
        return;
    }
    mpz_class d = pollard_rho(n);
    factor_recursive(d, factors);
    factor_recursive(n / d, factors);
}

static std::vector<mpz_class> distinct_prime_factors(
    const mpz_class &n) {
    mpz_class m = absz(n);
    if (m <= 1) return {};
    std::vector<mpz_class> factors;
    factor_recursive(m, factors);
    std::sort(factors.begin(), factors.end());
    factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
    return factors;
}

// ============================================================
// Factor support for this curve family
//
// B = 32(n+3)
//
// A² - 4B = (2n-3)(2n+5)^3
//
// Hence the relevant prime support is
//
//   2 * (n+3) * (2n-3) * (2n+5)
// ============================================================

static std::vector<mpz_class> descent_primes(const mpz_class &n) {
    std::vector<mpz_class> result;
    result.push_back(2);
    auto append = [&](const mpz_class &x) {
        auto f = distinct_prime_factors(x);
        result.insert(result.end(), f.begin(), f.end());
    };
    append(n + 3);
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

// ============================================================
// PARI rational parsing
// ============================================================

static mpq_class parse_rational(const std::string &s) {
    std::size_t slash = s.find('/');
    if (slash == std::string::npos)
        return mpq_class(mpz_class(s));
    mpz_class num(s.substr(0, slash));
    mpz_class den(s.substr(slash + 1));
    return mpq_class(num, den);
}

static std::string rational_to_string(const mpq_class &q) {
    std::ostringstream out;
    out << q.get_num();
    if (q.get_den() != 1)
        out << '/' << q.get_den();
    return out.str();
}
static Point dummy({0, -1});

static std::vector<Point> parse_pari_points(const std::string &input) {
    std::vector<Point> result;
    if (input.substr(0, 4) == "[[]]") return {dummy};
    std::size_t pos = input.find('[');
    while (pos != std::string::npos) {
        std::size_t open = input.find('[', pos + 1);
        if (open == std::string::npos) break;
        std::size_t close = input.find(']', open);
        if (close == std::string::npos) throw std::runtime_error("missing ]");
        std::string inside = input.substr(open + 1, close - open - 1);
        std::size_t comma = inside.find(',');
        if (comma != std::string::npos) {
            auto trim = [](std::string s) {
                const char *ws = " \t\r\n";
                std::size_t first = s.find_first_not_of(ws);
                std::size_t last = s.find_last_not_of(ws);
                if (first == std::string::npos) return std::string();
                return s.substr(first, last - first + 1);
            };
            std::string xs = trim(inside.substr(0, comma));
            std::string ys = trim(inside.substr(comma + 1));
            result.push_back({parse_rational(xs), parse_rational(ys)});
        }
        pos = close;
    }
    return result;
}

// ============================================================
// Construct d from valuation parity
// ============================================================

/*
 * Given a numerator a and the primes dividing
 *
 *     B * (A² - 4B),
 *
 * construct the squarefree part of a.
 *
 * We repeatedly divide by p and count v_p(a).
 *
 * If v_p(a) is odd, p is included in d.
 *
 * The sign of a is included separately.
 */
static mpz_class descent_squareclass(const mpz_class &a, const std::vector<mpz_class> &primes, mpz_class &square_part) {
    if (a == 0) {
        square_part = 0;
        return 0;
    }
    mpz_class m = absz(a);
    mpz_class d = a < 0 ? -1 : 1;
    for (const mpz_class &p : primes) {
        unsigned valuation = 0;
        while (m % p == 0) {
            m /= p;
            ++valuation;
        }
        if (valuation & 1) d *= p;
    }
    mpz_div(m.get_mpz_t(), a.get_mpz_t(), d.get_mpz_t());
    mpz_sqrt(square_part.get_mpz_t(), m.get_mpz_t());
    if (square_part * square_part != m)
        throw std::runtime_error("numerator contains an odd valuation outside the descent prime support");
    return d;
}

// ============================================================
// Compress one point
// ============================================================

static CompressedPoint compress_point(const Point &P, const std::vector<mpz_class> &primes) {
    mpz_class a = P.x.get_num();
    mpz_class den = P.x.get_den();
    if (den < 0) {
        den = -den;
        a = -a;
    }
    mpz_class s;
    mpz_sqrt(s.get_mpz_t(), den.get_mpz_t());
    if (s * s != den)
        throw std::runtime_error("x denominator is not a square");
    mpz_class r;
    mpz_class d = descent_squareclass(a, primes, r);
    if (d * r * r != a) {
        throw std::runtime_error("internal squareclass error");
    }
    return {d, r, s};
}

// ============================================================
// Decompression
// ============================================================

static mpq_class decompress_x(const CompressedPoint &p) {
    return mpq_class(p.d * p.r * p.r, p.s * p.s);
}

/*
 * For x = d*r²/s²:
 *
 * y² =
 *
 * x³ + A x² + Bx
 *
 * =
 *
 * d*r² *
 * (d²*r⁴ + A*d*r²*s² + B*s⁴)
 * / s⁶.
 *
 * Thus
 *
 * y =
 * r * sqrt(
 *     d³*r⁴
 *   + A*d²*r²*s²
 *   + B*d*s⁴
 * ) / s³.
 *
 * Equivalently we can simply construct x and use
 * rational arithmetic, but the expression below
 * avoids enormous rational intermediates.
 */
static mpq_class recover_y(const Curve &C, const CompressedPoint &p) {
    const mpz_class &d = p.d;
    const mpz_class &r = p.r;
    const mpz_class &s = p.s;

    mpz_class r2 = r * r;
    mpz_class s2 = s * s;
    mpz_class s4 = s2 * s2;
    mpz_class r4 = r2 * r2;

    /*
     * N satisfies
     *
     * y² = r²*N / s⁶.
     */
    mpz_class N = d * d * d * r4 + C.A * d * d * r2 * s2 + C.B * d * s4;
    if (N < 0)
        throw std::runtime_error("negative y²");
    if (!mpz_perfect_square_p(N.get_mpz_t()))
        throw std::runtime_error("y is not rational");
    mpz_class root;
    mpz_sqrt(root.get_mpz_t(), N.get_mpz_t());
    mpz_class y_num = r * root;
    mpz_class y_den = s * s * s;
    return mpq_class(y_num, y_den);
}

// ============================================================
// Complete compression/decompression
// ============================================================

static std::string compress_pari(const std::string &input) {
    std::istringstream in(input);
    std::string n_string;
    in >> n_string;
    if (n_string.empty())
        throw std::runtime_error("missing n");
    mpz_class n(n_string);
    std::size_t bracket = input.find('[');
    if (bracket == std::string::npos)
        throw std::runtime_error("missing points");
    auto points = parse_pari_points(input.substr(bracket));
    /*
     * Factor the support once for all points.
     */
    auto primes = descent_primes(n);
    std::ostringstream out;
    out << std::hex << n << ' ' << points.size();
    if (points.size() == 1 && points[0] == dummy) {
        out << ' ' << 0;
    }
    else {
        for (const auto &point : points) {
            CompressedPoint p = compress_point(point, primes);
            out << ' ' << p.d << ' ' << p.r << ' ' << p.s;
        }
    }
    return out.str();
}
static std::string decompress_pari(const std::string &input) {
    std::istringstream in(input);
    std::string n_string;
    std::size_t count;
    in >> std::hex >> n_string >> count;
    if (!in)
        throw std::runtime_error("invalid compressed input");
    Curve curve(mpz_class{n_string, 16});
    std::ostringstream out;
    out << curve.n << " [";
    for (std::size_t i = 0; i < count; ++i) {
        CompressedPoint p;
        if (!(in >> p.d))
            throw std::runtime_error("truncated input");
        if (p.d == 0) {
            out << "[]]";
            return out.str();
        }
        if (!(in >> p.r >> p.s))
            throw std::runtime_error("truncated input");
        mpq_class x = decompress_x(p);
        mpq_class y = recover_y(curve, p);
        if (i) out << ", ";
        out << '[' << rational_to_string(x) << ", " << rational_to_string(y) << ']';
    }
    out << ']';
    return out.str();
}

// ============================================================
// Main
// ============================================================

int main(int argc, char **argv) {
    try {
        bool compress = true, batch = false;
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "-d") == 0) compress = false;
            if (strcmp(argv[i], "--batch") == 0) batch = true;
            if (strcmp(argv[i], "--single") == 0) batch = false;
        }
        if (!batch) {
            std::string input;
            std::getline(std::cin, input);
            if (compress) {
                std::cout << compress_pari(input) << '\n';
            }
            else {
                std::cout << decompress_pari(input) << '\n';
            }
        }
        else {
            std::string input;
            while (std::getline(std::cin, input)) {
                if (compress) {
                    std::cout << compress_pari(input) << '\n';
                }
                else {
                    std::cout << decompress_pari(input) << '\n';
                }
            }
        }
    }
    catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}