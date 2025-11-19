path = r'src/app/main_imgui.cpp'
with open(path, 'r', encoding='utf-8') as f:
    lines = f.read().splitlines()

stack = []
for i, l in enumerate(lines, start=1):
    s = l.strip()
    if s.startswith('#if') and not s.startswith('#ifndef') and not s.startswith('#ifdef'):
        stack.append((s, i))
    elif s.startswith('#ifdef') or s.startswith('#ifndef'):
        stack.append((s, i))
    elif s.startswith('#endif'):
        if stack:
            stack.pop()
        else:
            print('Unmatched #endif at line', i)

if stack:
    print('Unmatched #if-like directives:')
    for s, i in stack:
        print(' ', i, s)
else:
    print('All #if/#endif balanced')
