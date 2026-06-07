"""Sphinx helpers for EASYclimate gallery references."""

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.statemachine import StringList

from gallery_examples import GALLERY_EXAMPLES


class EasyClimateMiniGallery(Directive):
    """Expand stable example aliases before delegating to sphinx-gallery."""

    has_content = True
    option_spec = {
        "add-heading": directives.unchanged,
        "heading-level": directives.nonnegative_int,
    }

    def run(self):
        gallery_examples = self.state.document.settings.env.config.easyclimate_gallery_examples
        source = self.state_machine.get_source(self.lineno - 1)
        lines = [".. minigallery::"]

        for name, value in self.options.items():
            if value is None:
                lines.append(f"   :{name}:")
            else:
                lines.append(f"   :{name}: {value}")

        lines.append("")
        for item in self.content:
            stripped = item.strip()
            resolved = gallery_examples.get(stripped, stripped)
            lines.append(f"   {resolved}" if stripped else "")

        node = nodes.container()
        self.state.nested_parse(StringList(lines, source=source), self.content_offset, node)
        return node.children


def setup(app):
    app.add_config_value("easyclimate_gallery_examples", GALLERY_EXAMPLES, "env")
    app.add_directive("ecl-minigallery", EasyClimateMiniGallery)
    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
