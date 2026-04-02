import csv
from itertools import islice

class Utils:
    @staticmethod
    def to_csv(iterator, filename, chunk_size=5000):
        try:
            first_row = next(iterator)
        except StopIteration:
            return

        with open(filename, 'w', newline='', buffering=128*1024, encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=first_row.keys())
            writer.writeheader()
            writer.writerow(first_row)
            while True:
                batch = list(islice(iterator, chunk_size))
                if not batch:
                    break
                writer.writerows(batch)

