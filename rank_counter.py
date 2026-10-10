from collections import Counter
from fractions import Fraction
import re, sys
from sage.all import magma
magma.eval('''
CheckLocalSolubility := function(n)
    P<x> := PolynomialRing(Rationals());
    A := 4*n^2 + 12*n - 3;
    B := 32*n + 96;
    for d in Divisors(SquareFreeFactorization(B)) do
        HyperE := HyperellipticCurve(-d*x^4 + A*x^2 - B/d);
        if IsLocallySoluble(GenusOneModel(HyperE)) and #FourDescent(HyperE) gt 0 then
            return true;
        end if;
    end for;
    return false;
end function;
''')
def check_curve(n):
    return bool(magma.CheckLocalSolubility(n))
filename = input("Database file: ").strip()
N = int(input("Number of lines: "))
def extract_bracket_block(s):
    start = s.find('[')
    if start == -1:
        return None
    depth = 0
    for i in range(start, len(s)):
        if s[i] == '[':
            depth += 1
        elif s[i] == ']':
            depth -= 1
            if depth == 0:
                return s[start:i + 1]
    return None
def get_rank(block):
    if block == "[]":
        return 0
    inner = block[1:-1].strip()
    if not inner:
        return 0
    return inner.count('[')
def get_x_values(block):
    if block == "[]":
        return []
    pairs = re.findall(
        r'\[\s*([^,\[\]]+)\s*,\s*([^,\[\]]+)\s*\]',
        block
    )
    xs = []
    for x, y in pairs:
        try:
            xs.append(Fraction(x))
        except ValueError:
            pass
    return xs
rank_count = Counter()
rank1_definite_solution = 0
rank1_potential_solution = 0
rank1_definite_no_solution = 0
rank_ge2_definite_solution = Counter()
generator_not_found = 0
invalid = 0
processed = 0
def print_aggregate_data():
    print(f"Processed: {processed}")
    print(f"Number of curves leading to a solution: "\
        f"{rank1_definite_solution + sum(rank_ge2_definite_solution[r] for r in sorted(rank_count)[2:])} + [0,{rank1_potential_solution}]")
    print(f"Generator not found: {generator_not_found}")
    print(f"Invalid: {invalid}")
    for r in sorted(rank_count):
        print(f"Rank {r}: {rank_count[r]} curves")
        if r == 1:
            print(
                f"- {rank1_definite_solution} curves "
                f"definitely lead to one solution"
            )
            print(
                f"- {rank1_potential_solution} curves "
                f"potentially lead to one solution"
            )
            print(
                f"- {rank1_definite_no_solution} curves "
                f"definitely don't lead to one solution"
            )
        elif r >= 2:
            print(
                f"({rank_ge2_definite_solution[r]} curves "
                f"definitely lead to one solution)"
            )
    print()
with open(filename, "r") as f:
    for _ in range(N):
        line = f.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        m = re.match(r'^\s*(\d+)(?:\s+(.*))?$', line)
        if not m:
            print(f'Invalid line: {line}')
            invalid += 1
            continue
        n = int(m.group(1))
        # Important:
        # If there is nothing after n, this is:
        #
        #     n
        #
        # which means Rank 1, generator not found.
        rest = m.group(2)
        processed += 1
        if rest is None or not rest.strip():
            rest = "[[]]"
        rest = rest.strip()
        block = extract_bracket_block(rest)
        if block is None:
            print(f'Invalid line: {line}')
            invalid += 1
            continue
        rank = get_rank(block)
        rank_count[rank] += 1
        if rank == 1:
            if block == "[[]]":
                generator_not_found += 1
                if n % 2 == 0 and check_curve(n):
                    rank1_potential_solution += 1
                else:
                    rank1_definite_no_solution += 1
                continue
            xs = get_x_values(block)
            if any(x < 0 for x in xs):
                rank1_definite_solution += 1
            else:
                rank1_definite_no_solution += 1
        elif rank >= 2:
            xs = get_x_values(block)
            if any(x < 0 for x in xs):
                rank_ge2_definite_solution[rank] += 1
        if n % 5000 == 0:
            print_aggregate_data()
