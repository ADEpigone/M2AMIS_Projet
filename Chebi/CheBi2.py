import sqlite3

class CheBi2:
    def __init__(self, bd_link = None):
        self.bd_link = bd_link
        self.conn = None
        self.cursor = None
        if self.bd_link:
            self.conn = sqlite3.connect(self.bd_link)
            self.cursor = self.conn.cursor()
            self._create_tables()
        print(f"Connecté à {self.bd_link}")

    def _create_tables(self):
        if not self.cursor:
            raise Exception("Database connection is not established.")

        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS chebi2 (
                chebi_id TEXT PRIMARY KEY,
                mol_name TEXT,
                mol_file TEXT
            )
        ''')
        self.conn.commit()

    def update_table(self, batch):
        if not self.cursor:
            raise Exception("Database connection is not established.")

        self.cursor.executemany('INSERT OR REPLACE INTO chebi2 (chebi_id, mol_name, mol_file) VALUES (?, ?, ?)', batch)
        self.conn.commit()

    def get_mol(self, chebi_id):
        if not self.cursor:
            raise Exception("Database connection is not established.")

        self.cursor.execute('SELECT mol_file FROM chebi2 WHERE chebi_id = ?', (chebi_id,))
        result = self.cursor.fetchone()
        if result:
            return result[0]
        return None

if __name__ == "__main__":
    bd = CheBi2("chebi2.db")
    bd.update_table([("CHEBI:37", "(+)-Nornicotine", "molfile content here")])
    print(bd.get_mol("CHEBI:37"))

