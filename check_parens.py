import io
path = r'src/app/main_imgui.cpp'
with open(path, 'r', encoding='utf-8') as f:
    text = f.read()

state = 'code'
line = 1
col = 0
stack = []  # list of (char, line, col)

i = 0
while i < len(text):
    ch = text[i]
    if ch == '\n':
        line += 1
        col = 0
        i += 1
        continue
    col += 1
    nxt = text[i+1] if i+1 < len(text) else ''

    if state == 'code':
        if ch == '/' and nxt == '/':
            state = 'line_comment'
            i += 2
            col += 1
            continue
        if ch == '/' and nxt == '*':
            state = 'block_comment'
            i += 2
            col += 1
            continue
        if ch == '"':
            state = 'string'
            i += 1
            continue
        if ch == '\'':
            state = 'char'
            i += 1
            continue
        if ch == '(':
            stack.append(('(', line, col))
        elif ch == ')':
            if stack and stack[-1][0] == '(':
                stack.pop()
            else:
                print('Unmatched ) at', line, col)
        i += 1
        continue

    if state == 'string':
        if ch == '\\':
            i += 2
            col += 1
            continue
        if ch == '"':
            state = 'code'
        i += 1
        continue

    if state == 'char':
        if ch == '\\':
            i += 2
            col += 1
            continue
        if ch == '\'':
            state = 'code'
        i += 1
        continue

    if state == 'line_comment':
        if ch == '\n':
            state = 'code'
        i += 1
        continue

    if state == 'block_comment':
        if ch == '*' and nxt == '/':
            state = 'code'
            i += 2
            col += 1
            continue
        i += 1
        continue

if stack:
    print('Unmatched ( positions:')
    for _, ln, cl in stack:
        print('  at line', ln, 'col', cl)
else:
    print('All parentheses matched')
