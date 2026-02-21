import argparse

from Chebi.CheBi import Chebi
from Chebi.CheBi2 import CheBi2
from utils import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    commands = parser.add_subparsers(dest='command', required=True)

    chebi_client = CheBi2("chebi2.db")

    plugins = [plugin(commands, chebi_client=chebi_client) for plugin in get_all_plugins()]


    args = parser.parse_args()

    for plugin in plugins:
        if plugin.match(args.command):
            plugin.execute(args)

    print(args)
