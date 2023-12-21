import random
import string
import pickle

# Funktion, um zufällige Zeichenketten zu generieren
def generate_random_string(length):
    letters = string.ascii_letters
    return ''.join(random.choice(letters) for _ in range(length))

# Anzahl der Strings und Länge jedes Strings
num_strings = 1
string_length = 4000000

# Liste mit zufälligen Zeichenketten generieren
string_list = [generate_random_string(string_length) for _ in range(num_strings)]

# Zeichenketten in Datei mit pickle speichern
with open("string_list.txt", "wb") as file:
    pickle.dump(string_list, file)

print(f"Die Zeichenketten wurden in string_list.txt gespeichert.")
