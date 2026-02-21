from cli_plugins.base.CLI_plugin import CLIPlugin
from analysis.correlation import run_correlation_benchmark


class CorrelationPlugin(CLIPlugin):

	def __init__(self, parser, **kwargs):
		super().__init__(
			parser,
			command_name="correlation",
			help_text="Lance un benchmark de corrélation syntaxique/sémantique",
		)
		self.add_argument(
			"--per-family",
			help_text="Nombre max de molécules par famille",
			default=40,
			type=int,
		)
		self.add_argument(
			"--kernel-similarity",
			help_text="Similarité du kernel CWL (tanimoto, dice, cosine)",
			default="cosine",
			choices=["tanimoto", "dice", "cosine"],
		)
		self.add_argument(
			"--fingerprint-method",
			help_text="Méthode de fingerprint syntaxique (cwl, morgan, rdkit)",
			default="cwl",
			choices=["cwl", "morgan", "rdkit"],
		)

	def execute(self, namespace, **kwargs):
		run_correlation_benchmark(
			molecules_per_class=namespace.per_family,
			kernel_similarity=namespace.kernel_similarity,
			fingerprint_method=namespace.fingerprint_method,
		)

