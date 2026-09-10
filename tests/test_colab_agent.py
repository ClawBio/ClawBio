import importlib.util
import json
import os
import tempfile
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('colab_agent', Path(__file__).resolve().parents[1] / 'examples/colab_agent.py')
agent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agent)

class AgentTests(unittest.TestCase):
    def test_real_dispatch_and_feedback(self):
        replies = iter([{'tool':'read_skill','args':{}}, {'tool':'run_analysis','args':{'sample':'baseline'}}, {'final':'Done'}])
        prompts, calls, events = [], [], []
        def model(prompt):
            prompts.append(prompt)
            return json.dumps(next(replies))
        session = agent.Session(model, {'read_skill':lambda: 'SKILL CONTENT', 'run_analysis':lambda sample: calls.append(sample) or {'value':42}}, events.append)
        session.ask('Analyse demo')
        self.assertEqual(calls, ['baseline'])
        self.assertIn('SKILL CONTENT', prompts[1])
        self.assertIn('42', prompts[2])
        self.assertEqual(events[-1]['kind'], 'answer')

    def test_unknown_tool_and_bad_json_recover(self):
        replies = iter(['not json', '{"tool":"shell","args":{"cmd":"bad"}}', '{"final":"Cannot do that"}'])
        session = agent.Session(lambda p:next(replies), {}, lambda e:None)
        session.ask('test')
        self.assertIn('Unknown tool', json.dumps(session.history))

    def test_limit_and_skill_gate(self):
        calls=[]
        session = agent.Session(lambda p:'{"tool":"run_analysis","args":{"sample":"baseline"}}', {'run_analysis':lambda sample:calls.append(sample)}, lambda e:None, max_steps=2)
        session.ask('run')
        self.assertEqual(calls, [])
        self.assertIn('step limit', json.dumps(session.history))

class NotebookTests(unittest.TestCase):
    def test_embedded_loop_matches_source(self):
        root = Path(__file__).resolve().parents[1]
        notebook = json.loads((root / 'docs/tutorial-agent.ipynb').read_text())
        code = [''.join(c['source']) for c in notebook['cells'] if c['cell_type'] == 'code']
        self.assertIn((root / 'examples/colab_agent.py').read_text().strip(), '\n'.join(code))
        for source in code:
            compile(source, 'notebook', 'exec')

if __name__ == '__main__': unittest.main()


class LaunchExperienceTests(unittest.TestCase):
    def test_single_hidden_launch_cell(self):
        root = Path(__file__).resolve().parents[1]
        notebook = json.loads((root / 'docs/tutorial-agent.ipynb').read_text())
        cells = [c for c in notebook['cells'] if c['cell_type'] == 'code']
        self.assertEqual(len(cells), 1)
        self.assertEqual(cells[0]['metadata']['cellView'], 'form')
        self.assertTrue(''.join(cells[0]['source']).startswith('#@title Launch ClawBio'))
