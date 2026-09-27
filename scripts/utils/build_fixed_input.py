#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用修复后的逻辑构建小输入，跑 ViReMa 验证不再 KeyError"""
import gzip

# 模拟修复后的 prep：取 R1/R2 各前 10000 条，tag 插入 name token
def build(out, limit=10000):
    n = 0
    with open(out, 'w') as o:
        for src, tag in [('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_1.fa.gz', '/1'),
                         ('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_2.fa.gz', '/2')]:
            cnt = 0
            with gzip.open(src, 'rt') as f:
                for line in f:
                    if cnt >= limit: break
                    line = line.rstrip('\n\r')
                    if line.startswith('>'):
                        parts = line.split(None, 1)
                        name = parts[0]
                        rest = (' ' + parts[1]) if len(parts) > 1 else ''
                        if not name.endswith(tag):
                            name = name + tag
                        o.write(name + rest + '\n')
                        cnt += 1
            n += cnt
    return n

n = build('/tmp/virema_test/fixed_input.fa')
print('reads written:', n)
print(open('/tmp/virema_test/fixed_input.fa').readline().rstrip())
