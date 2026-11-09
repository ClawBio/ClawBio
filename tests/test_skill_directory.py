"""Checks for the generated public directory."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('directory_generator', ROOT / 'scripts/generate_catalog.py')
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


def test_missing_directory_marker_fails():
    with pytest.raises(ValueError, match='exactly one'):
        generator.replace_directory_entries('<div>Changed template</div>', '<article>New</article>')


def test_directory_replacement():
    assert generator.replace_directory_entries(
        '<div id="results">old</div><!-- skills:end -->', 'new'
    ) == '<div id="results">new</div><!-- skills:end -->'


def test_every_skill_is_in_static_html():
    from html import escape
    skills = json.loads((ROOT / 'skills/directory.json').read_text())['skills']
    page = (ROOT / 'skills.html').read_text()
    assert page.count('<article ') == len(skills)
    for skill in skills:
        assert 'id="' + escape(skill['slug'], quote=True) + '"' in page
        assert '<h2 class="name">' + escape(skill['name']) + '</h2>' in page
        assert escape(skill['description']) in page
    assert 'No keywords specified' not in page
    assert '<span class="keywords"></span>' not in page


def test_only_tags_are_visible():
    from html import escape
    from html.parser import HTMLParser
    class Tags(HTMLParser):
        def __init__(self):
            super().__init__()
            self.skill = None
            self.in_tag = False
            self.tags = {}
        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if tag == 'article':
                self.skill = attrs['id']
                self.tags[self.skill] = []
            if tag == 'span' and attrs.get('class') == 'keyword':
                self.in_tag = True
        def handle_data(self, data):
            if self.in_tag:
                self.tags[self.skill].append(data)
        def handle_endtag(self, tag):
            if tag == 'span':
                self.in_tag = False
    parser = Tags()
    parser.feed((ROOT / 'skills.html').read_text())
    for skill in json.loads((ROOT / 'skills/directory.json').read_text())['skills']:
        assert parser.tags[skill['slug']] == list(dict.fromkeys(skill['tags']))
