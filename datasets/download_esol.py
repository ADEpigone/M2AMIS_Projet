import csv
import os
import urllib.request

DEFAULT_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv"


def _select_log_s_key(fieldnames):
    candidates = [
        "logS",
        "log_s",
        "measured log solubility in mols per litre",
    ]
    for key in candidates:
        if key in fieldnames:
            return key
    return None


def download_esol(output_path=None, url=DEFAULT_URL):
    if output_path is None:
        output_path = os.path.join("datasets", "esol.csv")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tmp_path = output_path + ".tmp"

    print(f"Telechargement ESOL depuis {url}")
    urllib.request.urlretrieve(url, tmp_path)

    with open(tmp_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise RuntimeError("CSV ESOL sans en-tete.")
        log_s_key = _select_log_s_key(reader.fieldnames)
        if log_s_key is None:
            raise RuntimeError("Colonne logS introuvable dans le CSV ESOL.")

        with open(output_path, "w", encoding="utf-8", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=["name", "smiles", "logS"])
            writer.writeheader()
            for row in reader:
                smiles = row.get("smiles") or row.get("SMILES")
                log_s = row.get(log_s_key)
                name = row.get("name") or row.get("Compound ID") or ""
                if smiles is None or log_s is None:
                    continue
                writer.writerow({"name": name, "smiles": smiles, "logS": log_s})

    os.remove(tmp_path)
    print(f"ESOL sauvegarde dans {output_path}")


if __name__ == "__main__":
    download_esol()
