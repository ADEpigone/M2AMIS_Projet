from datetime import timedelta
import sqlite3
import sys
sys.path.append("../")
from utils import get_mol_file
from datetime import datetime
class Chebi:
    def __init__(self, bd_link = "chebi_cache.db"):
        #sqlite bd
        self.bd_link = bd_link
        self.conn = None
        self.cursor = None
        if self.bd_link:
            self.conn = sqlite3.connect(self.bd_link)
            self.cursor = self.conn.cursor()
            self._create_tables()
        
    def _create_tables(self):
        if not self.cursor:
            raise Exception("Database connection is not established.")
        
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS chebi (
                chebi_id TEXT PRIMARY KEY,
                mol_file TEXT,
                date_added TEXT
            )
        ''')
        self.conn.commit()
    
    def get_mol(self, chebi_id):
        if not self.cursor:
            raise Exception("Database connection is not established.")
        
        self.cursor.execute('SELECT mol_file, date_added FROM chebi WHERE chebi_id = ?', (chebi_id,))
        result = self.cursor.fetchone()
        if result:
            mol = result[0]
            date = result[1]
            #si la date date d'il y a moins de 5j ? (arbitraire) on retourne
            #sinon on fait la récup comme après
            if datetime.now() - datetime.strptime(date, "%Y-%m-%d") < timedelta(days=5):
                return mol
               
        mol = get_mol_file(chebi_id)
        if mol:
            #on l'insère
            self.cursor.execute('REPLACE INTO chebi (chebi_id, mol_file, date_added) VALUES (?, ?, DATE("now"))',
                                (chebi_id, mol))
            self.conn.commit()

        return mol


if __name__ == "__main__":
    chebi = Chebi("chebi_cache.db")
    mol = chebi.get_mol("15377")
    print(mol)
    mol = chebi.get_mol("15377")




    