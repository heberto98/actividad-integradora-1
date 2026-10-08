# Actividad Integradora 1 - Parte 3: proteinas
# Traducimos el genoma en sus 3 marcos de lectura y buscamos cada proteina con KMP.
# La nsp12 tiene un frameshift (-1) a la mitad, asi que si una proteina no aparece
# completa la buscamos partida en dos marcos.

GENOMA = "SARS-COV-2-MN908947.3.txt"
PROTEINAS = "seq-proteins.txt"
BASE = 1  # indices en base 1 como en NCBI

NOMBRES = {
    "QHD43415_1": "nsp1", "QHD43415_2": "nsp2", "QHD43415_3": "nsp3",
    "QHD43415_4": "nsp4", "QHD43415_5": "nsp5 (3CLpro)", "QHD43415_6": "nsp6",
    "QHD43415_7": "nsp7", "QHD43415_8": "nsp8", "QHD43415_9": "nsp9",
    "QHD43415_10": "nsp10", "QHD43415_11": "nsp12 (RdRp)", "QHD43415_12": "nsp13 (helicasa)",
    "QHD43415_13": "nsp14", "QHD43415_14": "nsp15", "QHD43415_15": "nsp16",
    "QHD43416": "S (Spike)", "QHD43417": "ORF3a", "QHD43418": "E (envoltura)",
    "QHD43419": "M (membrana)", "QHD43420": "ORF6", "QHD43421": "ORF7a",
    "QHD43422": "ORF8", "QHD43423": "N (nucleocapside)", "QHI42199": "ORF10",
}

# codigo genetico estandar, * = stop
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codones = {}
n = 0
for a in bases:
    for b in bases:
        for c in bases:
            codones[a + b + c] = aas[n]
            n += 1


def leer_genoma(ruta):
    sec = []
    with open(ruta) as f:
        for linea in f:
            if linea.startswith(">"):
                continue
            sec.append(linea.strip().upper())
    return "".join(sec)


def leer_proteinas(ruta):
    lista = []
    nombre = None
    sec = []
    with open(ruta) as f:
        for linea in f:
            linea = linea.strip()
            if linea == "":
                continue
            if linea.startswith(">"):
                if nombre is not None:
                    lista.append((nombre, "".join(sec)))
                nombre = linea[1:].split()[0]
                sec = []
            else:
                sec.append(linea)
    if nombre is not None:
        lista.append((nombre, "".join(sec)))
    return lista


def traducir(genoma, marco):
    # el aminoacido k sale del codon que empieza en marco + 3k
    res = []
    for i in range(marco, len(genoma) - 2, 3):
        res.append(codones.get(genoma[i:i + 3], "X"))
    return "".join(res)


def invertir(s):
    res = []
    for i in range(len(s) - 1, -1, -1):
        res.append(s[i])
    return "".join(res)


# ---------- KMP ----------

def lps(p):
    tabla = [0] * len(p)
    k = 0
    for i in range(1, len(p)):
        while k > 0 and p[i] != p[k]:
            k = tabla[k - 1]
        if p[i] == p[k]:
            k += 1
        tabla[i] = k
    return tabla


def kmp(texto, p, tabla):
    # regresa las posiciones donde aparece p en texto
    pos = []
    k = 0
    for i in range(len(texto)):
        while k > 0 and texto[i] != p[k]:
            k = tabla[k - 1]
        if texto[i] == p[k]:
            k += 1
        if k == len(p):
            pos.append(i - k + 1)
            k = tabla[k - 1]
    return pos


def kmp_prefijos(texto, p, tabla):
    # pref[i] = que tan largo es el prefijo de p que termina en texto[i]
    pref = [0] * len(texto)
    k = 0
    for i in range(len(texto)):
        while k > 0 and (k == len(p) or texto[i] != p[k]):
            k = tabla[k - 1]
        if texto[i] == p[k]:
            k += 1
        pref[i] = k
    return pref


# ---------- busqueda ----------
# Cada resultado se guarda como un diccionario:
#   m = largo de la proteina
#   los primeros k aminoacidos salen del marco1 desde la posicion ini1
#   los demas (si k < m) salen del marco2 desde la posicion ini2

def pos_codon(r, t):
    # posicion en el genoma del codon que produce el aminoacido t
    if t < r["k"]:
        return r["marco1"] + 3 * (r["ini1"] + t)
    return r["marco2"] + 3 * (r["ini2"] + t - r["k"])


def rangos(r):
    a = r["marco1"] + 3 * r["ini1"]
    res = [(a, a + 3 * r["k"] - 1)]
    if r["k"] < r["m"]:
        b = r["marco2"] + 3 * r["ini2"]
        res.append((b, b + 3 * (r["m"] - r["k"]) - 1))
    return res


def buscar(prot, marcos):
    tabla = lps(prot)
    m = len(prot)
    res = []
    for f in range(3):
        for p in kmp(marcos[f], prot, tabla):
            res.append({"m": m, "k": m, "marco1": f, "ini1": p, "marco2": f, "ini2": p + m})
    return res


def buscar_frameshift(prot, marcos):
    # Buscamos un pedazo inicial de la proteina en un marco y el pedazo final en otro.
    # Prefijos: KMP normal. Sufijos: KMP con la proteina y el texto al reves.
    m = len(prot)
    tabla = lps(prot)
    prot_inv = invertir(prot)
    tabla_inv = lps(prot_inv)

    prefijos = {}  # (marco, inicio) -> largo
    for f in range(3):
        t = marcos[f]
        pref = kmp_prefijos(t, prot, tabla)
        for i in range(len(t)):
            L = pref[i]
            if L > 0 and (i == len(t) - 1 or pref[i + 1] != L + 1):
                prefijos[(f, i - L + 1)] = L

    sufijos = {}  # (marco, fin) -> largo
    for f in range(3):
        t = invertir(marcos[f])
        pref = kmp_prefijos(t, prot_inv, tabla_inv)
        for i in range(len(t)):
            L = pref[i]
            if L > 0 and (i == len(t) - 1 or pref[i + 1] != L + 1):
                fin = len(t) - 1 - (i - L + 1)
                if L > sufijos.get((f, fin), 0):
                    sufijos[(f, fin)] = L

    res = []
    for (a, ini), L1 in prefijos.items():
        for b in range(3):
            # d = distancia entre el ultimo nucleotido del pedazo 1 y el primero del pedazo 2
            # d = 1 seria seguir normal, d = 0 es frameshift -1 y d = 2 es +1
            for d in (0, 2):
                x = d - 1 - b + a
                if x % 3 != 0:
                    continue
                fin = ini + m - 1 + x // 3
                L2 = sufijos.get((b, fin), 0)
                if L2 > 0 and L1 + L2 >= m:
                    res.append({"m": m, "k": L1, "marco1": a, "ini1": ini,
                                "marco2": b, "ini2": fin - (m - L1) + 1})
    return res


def parte3():
    print("===== PARTE 3: PROTEINAS =====")
    genoma = leer_genoma(GENOMA)
    proteinas = leer_proteinas(PROTEINAS)
    marcos = [traducir(genoma, 0), traducir(genoma, 1), traducir(genoma, 2)]

    encontradas = 0
    faltantes = []
    for nombre, prot in proteinas:
        res = buscar(prot, marcos)
        shift = False
        # en proteinas muy cortas un prefijo + sufijo podria salir por casualidad
        if not res and len(prot) >= 20:
            res = buscar_frameshift(prot, marcos)
            shift = len(res) > 0
        if not res:
            faltantes.append(nombre)
            continue

        encontradas += 1
        titulo = nombre
        if nombre in NOMBRES:
            titulo = nombre + " - " + NOMBRES[nombre]
        print(f"\n{titulo} ({len(prot)} aa)")
        for r in res:
            texto = []
            for a, b in rangos(r):
                texto.append(f"{a + BASE} - {b + BASE}")
            extra = "  (frameshift)" if shift else ""
            print("  Indices:", ", ".join(texto) + extra)
            cod = []
            for t in range(min(4, len(prot))):
                p = pos_codon(r, t)
                cod.append(genoma[p:p + 3])
            print("  Primeros 4 aa:", prot[:4])
            print("  Codones:", " ".join(cod))

    print(f"\nProteinas encontradas: {encontradas} de {len(proteinas)}")
    if faltantes:
        print("No encontradas:", ", ".join(faltantes))


if __name__ == "__main__":
    parte3()
