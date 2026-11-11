import argparse
import json
import sys
from pathlib import Path
from .records import inspect_bundle, load_rows
from .segmentation import process_bundle
from .storage import RevisionStore, confined, save_json


def run_segment(bundle_root, output_path):
    """读取资料目录中的三类记录并执行离线空间归一与裁切。"""
    root = Path(bundle_root)
    transforms = load_rows(confined(root, 'transforms.json'))
    regions = load_rows(confined(root, 'regions.json'))
    paths = load_rows(confined(root, 'paths.json'))
    report = process_bundle(transforms, regions, paths)
    if output_path is not None:
        target = output_path if output_path.is_absolute() \
            else (Path.cwd() / output_path).resolve()
        save_json(target, report)
    return report


def main():
    parser = argparse.ArgumentParser(description='数据包内容检查与空间裁切')
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('inspect')
    p.add_argument('bundle', type=Path)
    p = commands.add_parser('catalog')
    p.add_argument('database', type=Path)
    p.add_argument('stream')
    p = commands.add_parser('segment', description='离线空间归一与边界裁切')
    p.add_argument('bundle', type=Path)
    p.add_argument('--output', type=Path, default=None, help='可选：报告写出路径')
    args = parser.parse_args()
    try:
        if args.command == 'inspect':
            value = inspect_bundle(args.bundle)
        elif args.command == 'catalog':
            value = RevisionStore(args.database).latest(args.stream)
        else:
            value = run_segment(args.bundle, args.output)
        print(json.dumps(value, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, KeyError, OSError):
        print('输入资料或目录无效，请检查路径及字段说明', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
