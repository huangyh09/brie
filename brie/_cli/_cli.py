import click

from .count import count
from .quant import quant
from .extract import extract
from .download import download

def _get_version():
    import pkg_resources
    import re
    from os.path import realpath, dirname, join

    if __name__ == "__main__":
        print('test1')
        filepath = join(dirname(realpath(__file__)), "..", "version.py")
        with open(filepath, "r") as f:
            content = f.read()
    else:
        content = pkg_resources.resource_string(__name__.split(".")[0], "version.py")
        content = content.decode()

    c = re.compile(r"__version__ *= *('[^']+'|\"[^\"]+\")")
    m = c.search(content)
    if m is None:
        return "unknown"
    return m.groups()[0][1:-1]


@click.group(name="brie", context_settings=dict(help_option_names=["-h", "--help"]))
@click.pass_context
@click.version_option(_get_version())
def cli(ctx):
    pass


cli.add_command(count)
cli.add_command(quant)
cli.add_command(extract)
cli.add_command(download)
