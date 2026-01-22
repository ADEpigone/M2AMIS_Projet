import argparse

from utils import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    commands = parser.add_subparsers(dest='command', required=True)

    plugins = [plugin(commands) for plugin in get_all_plugins()]


    args = parser.parse_args()

    for plugin in plugins:
        if plugin.match(args.command):
            plugin.execute(args)

    print(args)
