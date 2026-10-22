"""Keep every Python README example executable, in reading order."""

import re
from pathlib import Path


def test_readme_python_examples():
    readme = Path(__file__).parents[1] / "README.md"
    blocks = re.findall(r"```python\n(.*?)```", readme.read_text(), flags=re.DOTALL)
    assert blocks, "README must contain runnable Python examples"
    namespace = {}
    for index, source in enumerate(blocks):
        exec(compile(source, f"README.md:python-block-{index}", "exec"), namespace)
