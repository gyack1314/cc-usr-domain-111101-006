import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def test_real_inspection(self):
        result = subprocess.run([sys.executable, '-m', 'datapack.cli', 'inspect', 'samples'],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(json.loads(result.stdout)['files'])

    def test_invalid_path_is_failure(self):
        result = subprocess.run([sys.executable, '-m', 'datapack.cli', 'inspect',
                                 'missing-directory'], cwd=ROOT,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, '')

    def test_segment_report_on_samples(self):
        result = subprocess.run([sys.executable, '-m', 'datapack.cli', 'segment', 'samples'],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(report['summary']['total'], 7)
        self.assertEqual(report['summary']['accepted'], 5)
        self.assertEqual(report['summary']['rejected'], 2)
        accepted = [item for item in report['paths'] if item['accepted']]
        self.assertTrue(all(item['verification']['all_ok'] for item in accepted))
        # 缺变换链与冲突链各拒绝一条
        reasons = ' '.join(item['reason'] for item in report['paths']
                           if not item['accepted'])
        self.assertIn('无可达变换链', reasons)
        self.assertIn('不一致的变换链', reasons)

    def test_segment_writes_output_file(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'report.json'
            result = subprocess.run(
                [sys.executable, '-m', 'datapack.cli', 'segment', 'samples',
                 '--output', str(output)],
                cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(output.exists())
            report = json.loads(output.read_text(encoding='utf-8'))
            self.assertIn('summary', report)
