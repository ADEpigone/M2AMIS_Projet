import os

from cli_plugins.base.CLI_plugin import CLIPlugin
from analysis.comp_prediction import compare_prediction_benchmark
from analysis.comp_similarity import compare_similarity_kernels
from datasets.download_esol import download_esol
from utils import load_esol_dataset


class PredictionBenchmarkPlugin(CLIPlugin):

	def __init__(self, parser, **kwargs):
		super().__init__(
			parser,
			command_name="comp_prediction",
			help_text="Benchmark de prédiction logS avec différents fingerprints/modèles",
		)
		self.add_argument(
			"--operation",
			help_text="Opération à lancer: prediction, similarity, all",
			default="all",
			choices=["prediction", "similarity", "all"],
		)
		self.add_argument(
			"--dataset",
			help_text="Chemin du CSV ESOL",
			default=os.path.join("datasets", "esol.csv"),
		)
		self.add_argument(
			"--limit",
			help_text="Nombre max de molécules à charger (<=0 pour tout)",
			default=0,
			type=int,
		)
		self.add_argument(
			"--seed",
			help_text="Seed pour le mélange du dataset",
			default=42,
			type=int,
		)
		self.add_argument(
			"--fingerprints",
			help_text="Liste des fingerprints séparés par des virgules (morgan,rdkit,cwl)",
			default="morgan,rdkit,cwl",
		)
		self.add_argument(
			"--models",
			help_text="Liste des modèles séparés par des virgules (ridge,rf)",
			default="ridge,rf",
		)
		self.add_argument(
			"--folds",
			help_text="Nombre de folds de validation croisée",
			default=5,
			type=int,
		)
		self.add_argument(
			"--experiment-name",
			help_text="Nom de l'expérience",
			default="esol_prediction",
		)
		self.add_argument(
			"--similarities",
			help_text="Liste des similarités séparées par des virgules (tanimoto,dice,cosine)",
			default="tanimoto,dice,cosine",
		)
		self.add_argument(
			"--sample-frac",
			help_text="Fraction de données utilisée pour le benchmark similarité",
			default=0.8,
			type=float,
		)
		self.add_argument(
			"--target-prop",
			help_text="Propriété cible pour similarité",
			default="logS",
		)

	@staticmethod
	def _parse_csv_list(raw: str):
		values = [part.strip().lower() for part in str(raw).split(",") if part.strip()]
		return tuple(values)

	def execute(self, namespace, **kwargs):
		operation = namespace.operation
		dataset_path = namespace.dataset
		limit = None if namespace.limit <= 0 else namespace.limit
		seed = namespace.seed
		n_splits = namespace.folds

		fingerprint_names = self._parse_csv_list(namespace.fingerprints)
		model_names = self._parse_csv_list(namespace.models)
		similarity_names = self._parse_csv_list(namespace.similarities)
		similarity_builtin_fps = tuple(fp for fp in fingerprint_names if fp in {"morgan", "rdkit"})

		ignored_for_similarity = [fp for fp in fingerprint_names if fp not in {"morgan", "rdkit"}]
		if ignored_for_similarity and operation in {"similarity", "all"}:
			print(
				"Fingerprints ignorés pour similarity (non builtin): "
				+ ", ".join(ignored_for_similarity)
			)

		if not os.path.exists(dataset_path):
			print(f"Dataset introuvable: {dataset_path}")
			print("Téléchargement automatique d'ESOL...")
			download_esol(output_path=dataset_path)

		molecules, dataset_info = load_esol_dataset(dataset_path, limit=limit, seed=seed)
		print(f"Dataset ESOL: {dataset_info}")

		if not molecules:
			print("Dataset ESOL vide ou invalide. Impossible de lancer la prédiction.")
			return

		if operation in {"prediction", "all"} and n_splits < 2:
			print("--folds doit être >= 2")
			return

		if operation in {"prediction", "all"}:
			compare_prediction_benchmark(
				molecules,
				experiment_name=namespace.experiment_name,
				fingerprint_names=fingerprint_names,
				model_names=model_names,
				n_splits=n_splits,
			)

		if operation in {"similarity", "all"}:
			if not similarity_builtin_fps:
				print("Aucun fingerprint builtin valide pour similarity (utilise morgan et/ou rdkit).")
				return

			compare_similarity_kernels(
				molecules,
				prop=namespace.target_prop,
				experiment_name=f"{namespace.experiment_name}_similarity",
				similarity_names=similarity_names,
				builtin_fingerprints=similarity_builtin_fps,
				sample_frac=namespace.sample_frac,
			)

		print("Benchmark terminé.")

