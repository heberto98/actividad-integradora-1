# ============================================================
# E1. Actividad Integradora 1 - Parte 2
# Palíndromo más largo en cada gen (M, S y ORF1AB)
# Algoritmo: Manacher  ->  Tiempo O(n), Espacio O(n)
# ============================================================


def leer_secuencia(ruta):
    """Lee un archivo de secuencia. Ignora encabezados tipo FASTA ('>...')
    y conserva solo letras en mayúscula (quita saltos de línea, espacios, etc.)."""
    try:
        with open(ruta, "r") as archivo:
            letras = []
            for linea in archivo:
                if linea.startswith(">"):
                    continue
                for c in linea:
                    if "a" <= c <= "z":
                        c = chr(ord(c) - 32)
                    if "A" <= c <= "Z":
                        letras.append(c)
            return "".join(letras)
    except FileNotFoundError:
        print("Error: no se pudo abrir", ruta)
        return ""


def palindromo_mas_largo(s):
    """Algoritmo de Manacher.
    Se transforma s = "abc" en t = "^#a#b#c#$" para tratar igual los
    palíndromos de longitud par e impar. '^' y '$' son centinelas distintos
    que detienen la expansión sin revisar límites del arreglo.
    p[i] = radio del palíndromo centrado en t[i] = longitud del palíndromo en s.
    Regresa (inicio, longitud), con inicio en base 0."""
    n = len(s)
    if n == 0:
        return 0, 0

    m = 2 * n + 3
    t = ["#"] * m
    t[0] = "^"
    t[m - 1] = "$"
    for i in range(n):
        t[2 * i + 2] = s[i]

    p = [0] * m
    centro = 0
    derecha = 0          # límite derecho del palíndromo que llega más lejos
    mejor_centro = 0
    mejor_radio = 0

    for i in range(1, m - 1):
        # Aprovecha la simetría: el espejo de i respecto a 'centro'
        if i < derecha:
            espejo = p[2 * centro - i]
            p[i] = derecha - i if derecha - i < espejo else espejo

        # Expande mientras los caracteres coincidan
        while t[i + p[i] + 1] == t[i - p[i] - 1]:
            p[i] += 1

        # Actualiza el palíndromo que llega más a la derecha
        if i + p[i] > derecha:
            centro = i
            derecha = i + p[i]

        # Guarda el mejor
        if p[i] > mejor_radio:
            mejor_radio = p[i]
            mejor_centro = i

    # Convierte la posición en t a la posición en s
    inicio = (mejor_centro - mejor_radio - 1) // 2
    return inicio, mejor_radio


def parte2():
    genes = [
        ("M", "gen-M.txt"),
        ("S", "gen-S.txt"),
        ("ORF1AB", "gen-ORF1AB.txt"),
    ]

    print("===== PARTE 2: PALINDROMO MAS LARGO POR GEN =====")

    with open("palindromos.txt", "w") as salida:
        for nombre, ruta in genes:
            gen = leer_secuencia(ruta)
            if gen == "":
                continue

            inicio, longitud = palindromo_mas_largo(gen)
            texto = gen[inicio:inicio + longitud]

            # Posiciones mostradas en base 1 (como en NCBI)
            ini = inicio + 1
            fin = inicio + longitud

            print(f"\nGen {nombre} (longitud del gen: {len(gen)})")
            print(f"  Longitud del palindromo mas largo: {longitud}")
            print(f"  Posicion en el gen: {ini} - {fin}")
            print(f"  Palindromo: {texto}")

            salida.write(f"Gen {nombre}\n")
            salida.write(f"Longitud: {longitud}\n")
            salida.write(f"Posicion en el gen: {ini} - {fin}\n")
            salida.write(f"Palindromo: {texto}\n\n")

    print("\nPalindromos guardados en palindromos.txt")


if __name__ == "__main__":
    parte2()
