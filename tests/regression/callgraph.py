"""Static call graph of a Fortran program: units reachable from the main program.
    python3 callgraph.py MAIN.f90 file.f90 ...   -> reachable / unreachable units per file"""
import re, sys
files = sys.argv[1:]
units = {}      # name -> (file, body)
for f in files:
    txt = open(f, errors='replace').read()
    for m in re.finditer(r'^[ \t]*(?:(?:pure|elemental|recursive|integer|logical|real(?:\(\w+\))?|complex(?:\(\w+\))?|double precision|type\(\w+\))[ \t]+)*(subroutine|function|program)[ \t]+(\w+)(.*?)^[ \t]*end[ \t]*(?:\1)?[ \t]*\w*[ \t]*$', txt, re.M | re.S | re.I):
        units.setdefault(m.group(2).lower(), (f, m.group(3)))
names = set(units)
def callees(body):
    body = '\n'.join(l.split('!')[0] for l in body.splitlines())
    out = set(x.lower() for x in re.findall(r'\bcall\s+(\w+)', body, re.I))
    for n in names:                                  # function references
        if re.search(r'\b' + n + r'\s*\(', body, re.I): out.add(n)
    return out & names
main = [n for n, (f, b) in units.items() if f == files[0] and re.search(r'program', open(f).read()[:4000], re.I)]
start = [n for n in units if units[n][0] == files[0]]
seen, todo = set(), list(start)
while todo:
    n = todo.pop()
    if n in seen: continue
    seen.add(n); todo += list(callees(units[n][1]) - seen)
byfile = {}
for n, (f, b) in units.items():
    byfile.setdefault(f, [0, 0, []]); byfile[f][0 if n in seen else 1] += 1
    if n not in seen: byfile[f][2].append(n)
for f, (r, u, l) in sorted(byfile.items()):
    print(f'{f.split("/")[-1]:22s} reachable {r:4d}  unreachable {u:4d}')
print('total reachable', len(seen), 'of', len(units))
