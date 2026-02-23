import numpy as np

from Chebi.CheBi2 import CheBi2
from cli_plugins.base.CLI_plugin import CLIPlugin
from analysis.clusters import (
	DIST_THRESHOLD,
	JSON_OUTPUT,
	MAX_MOLECULES,
	load_clusters,
	run_clustering_and_save,
)
from analysis.dominante_families import (
	consensus_scores_for_all_clusters,
	plot_cumulative_curve,
)
from analysis.clustering_heatmap import generate_hierarchical_similarity_heatmap
from analysis.purity_analysis import analyze_family_purity
from similarites.builtin_similarity import BuiltinSimilarity
from similarites.cwl_kernel import CWLKernel
from similarites.ontology_similarity import OntologySimilarity
from utils import (
	MOL_LITE_UPDATE,
	MOL_LITE_URL,
	build_kernel,
	ensure_chebi_db_up_to_date,
	get_mol_lite,
	has_db_changed,
	load_ontology,
	prep_db_load,
)


class ClusteringPlugin(CLIPlugin):

	def __init__(self, parser, chebi_client: CheBi2, **kwargs):
		super().__init__(
			parser,
			command_name="clustering",
			help_text="Calcul et analyses de clustering sur la base ChEBI",
		)

		self.add_argument(
			"--operation",
			help_text="Opération à exécuter: run, load, purity, dominant, cdf, heatmap",
			required=True,
			choices=["run", "load", "purity", "dominant", "cdf", "heatmap"],
		)
		self.add_argument(
			"--input-file",
			help_text="Fichier JSON des clusters en entrée",
			default=JSON_OUTPUT,
		)
		self.add_argument(
			"--output-file",
			help_text="Fichier de sortie (run/cdf/heatmap)",
			default=None,
		)
		self.add_argument(
			"--kernel",
			help_text="Noyau de similarité pour run (cwl, morgan, rdkit, ontology)",
			default="cwl",
			choices=["cwl", "morgan", "rdkit", "ontology"],
		)
		self.add_argument(
			"--similarity",
			help_text="Mesure de similarité (tanimoto, dice, cosine)",
			default="tanimoto",
			choices=["tanimoto", "dice", "cosine"],
		)
		self.add_argument(
			"--threshold",
			help_text="Seuil de distance pour le clustering hiérarchique",
			default=DIST_THRESHOLD,
			type=float,
		)
		self.add_argument(
			"--sample-size",
			help_text="Taille de l'échantillon pour le clustering (<=0 pour tout)",
			default=MAX_MOLECULES,
			type=int,
		)
		self.add_argument(
			"--family",
			help_text="Famille ChEBI (nom ou CHEBI:ID) pour l'analyse de pureté",
			default=None,
		)
		self.add_argument(
			"--min-depth",
			help_text="Profondeur minimale ontologique",
			default=3,
			type=int,
		)
		self.add_argument(
			"--min-cluster-size",
			help_text="Taille minimale de cluster à analyser",
			default=2,
			type=int,
		)

		self.chebi_client = chebi_client

	def _load_clusters_or_exit(self, input_file: str):
		cluster_map = load_clusters(input_file)
		if cluster_map is None:
			return None
		return cluster_map


	def execute(self, namespace, **kwargs):
		operation = namespace.operation.lower()
		input_file = namespace.input_file
		output_file = namespace.output_file

		ensure_chebi_db_up_to_date()

		if operation == "run":
			sim_kernel, has_fingerprint = build_kernel(namespace.kernel, namespace.similarity)
			save_path = output_file if output_file else input_file
			max_molecules = None if namespace.sample_size <= 0 else namespace.sample_size
			run_clustering_and_save(
				sim_kernel=sim_kernel,
				save_path=save_path,
				has_fingerprint=has_fingerprint,
				dist_threshold=namespace.threshold,
				max_molecules=max_molecules,
			)
			return

		if operation == "load":
			cluster_map = self._load_clusters_or_exit(input_file)
			if cluster_map is None:
				return
			sizes = [len(members) for members in cluster_map.values()]
			print(f"Clusters chargés: {len(cluster_map)}")
			print(f"Taille moyenne: {np.mean(sizes):.2f} | Taille max: {max(sizes)}")
			return

		if operation == "purity":
			if not namespace.family:
				print("L'option --family est obligatoire pour l'opération purity.")
				return
			cluster_map = self._load_clusters_or_exit(input_file)
			if cluster_map is None:
				return
			ontology = load_ontology()
			molecules = list(self.chebi_client.get_all_mols())
			analyze_family_purity(namespace.family, ontology, molecules, cluster_map)
			return

		if operation == "dominant":
			cluster_map = self._load_clusters_or_exit(input_file)
			if cluster_map is None:
				return
			ontology = load_ontology()
			scores, details = consensus_scores_for_all_clusters(
				ontology,
				cluster_map,
				min_cluster_size=namespace.min_cluster_size,
				min_depth=namespace.min_depth,
			)
			if not scores:
				print("Aucun cluster exploitable pour l'analyse CDS des familles dominantes.")
				return
			score_arr = np.array(scores)
			print(f"Nombre de clusters analysés: {len(scores)}")
			print(f"Moyenne CDS : {score_arr.mean():.3f} | Médiane CDS : {np.median(score_arr):.3f}")
			print(f"% clusters avec CDS >= 0.5: {(np.mean(score_arr >= 0.5) * 100):.1f}%")

			top = sorted(
				details.items(),
				key=lambda item: item[1]["consensus_depth_score"],
				reverse=True)[:5]
			print("Top 5 clusters par score CDS:")
			for cid, info in top:
				print(
					f"  - Cluster {cid}: size={info['size']}, CDS={info['consensus_depth_score']:.3f}, "
					f"ratio={info['dominance_ratio']:.3f}, profondeur={info['depth']}, famille={info['dominant_family_id']}"
				)
			return

		if operation == "cdf":
			cluster_map = self._load_clusters_or_exit(input_file)
			if cluster_map is None:
				return
			ontology = load_ontology()
			scores, _ = consensus_scores_for_all_clusters(
				ontology,
				cluster_map,
				min_cluster_size=namespace.min_cluster_size,
				min_depth=namespace.min_depth,
			)
			if not scores:
				print("Aucun score CDS disponible pour tracer la courbe cumulative.")
				return

			plot_path = output_file if output_file else "cdf_consensus_depth_clusters.png"
			plot_cumulative_curve(
				scores,
				title="Courbe cumulative du score consensus-profondeur (CDS)",
				output_path=plot_path,
			)
			return

		if operation == "heatmap":
			sim_kernel, has_fingerprint = build_kernel(namespace.kernel, namespace.similarity)
			plot_path = output_file if output_file else "heatmap_hierarchical_similarity.png"
			max_molecules = None if namespace.sample_size <= 0 else namespace.sample_size
			generate_hierarchical_similarity_heatmap(
				sim_kernel=sim_kernel,
				has_fingerprint=has_fingerprint,
				dist_threshold=namespace.threshold,
				max_molecules=max_molecules,
				output_path=plot_path,
			)
			return

		print(f"Opération non supportée: {operation}")

