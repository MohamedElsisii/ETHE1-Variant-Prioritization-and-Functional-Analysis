
from pymol import cmd
from math import sqrt


# ============================================================
# SETTINGS
# ============================================================

# 9 ACTIVE-SITE RESIDUES
ACTIVE_SITE = [
    79,
    83,
    135,
    154,
    157,
    161,
    162,
    163,
    197
]

# 3 FE-BINDING RESIDUES
FE_BINDING = [
    79,
    135,
    154
]

# Maximum distance
MAX_DISTANCE = 6.0

# Distance used to distinguish stronger vs weaker association
STRONG_DISTANCE = 4.0

# Backbone atoms excluded from residue-residue contacts
BACKBONE = [
    "N",
    "CA",
    "C",
    "O"
]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_atoms(selection):

    return cmd.get_model(selection).atom


def distance(a, b):

    return sqrt(
        (a.coord[0] - b.coord[0]) ** 2 +
        (a.coord[1] - b.coord[1]) ** 2 +
        (a.coord[2] - b.coord[2]) ** 2
    )


def residue_name(resi):

    atoms = get_atoms(
        f"chain A and resi {resi}"
    )

    if len(atoms) > 0:
        return atoms[0].resn

    return "UNK"


def sidechain_atoms(resi):

    atoms = get_atoms(
        f"chain A and resi {resi}"
    )

    return [
        a for a in atoms
        if a.name not in BACKBONE
    ]


# ============================================================
# HEADER
# ============================================================

print("\n")
print("=" * 80)
print("ACTIVE-SITE + FE-BINDING STRUCTURAL NETWORK")
print("=" * 80)

print(
    f"Maximum distance = {MAX_DISTANCE:.1f} A"
)

print()

print("9 ACTIVE-SITE RESIDUES:")

print(
    " ".join(
        f"{residue_name(r)}{r}"
        for r in ACTIVE_SITE
    )
)

print()

print("3 FE-BINDING RESIDUES:")

print(
    " ".join(
        f"{residue_name(r)}{r}"
        for r in FE_BINDING
    )
)


# ============================================================
# GET ALL PROTEIN RESIDUES
# ============================================================

all_atoms = get_atoms(
    "chain A and polymer.protein"
)

all_residues = sorted(
    set(
        int(a.resi)
        for a in all_atoms
    )
)


# ============================================================
# ANALYZE ACTIVE-SITE CONNECTIONS
# ============================================================

active_network = {}


for reference_r in ACTIVE_SITE:

    reference_atoms = sidechain_atoms(
        reference_r
    )

    if not reference_atoms:
        continue


    for candidate_r in all_residues:

        # Don't analyze the reference residues themselves
        if candidate_r in ACTIVE_SITE:
            continue

        candidate_atoms = sidechain_atoms(
            candidate_r
        )

        if not candidate_atoms:
            continue


        best = None


        for ref_atom in reference_atoms:

            for candidate_atom in candidate_atoms:

                d = distance(
                    ref_atom,
                    candidate_atom
                )

                if d <= MAX_DISTANCE:

                    if (
                        best is None
                        or d < best[0]
                    ):

                        best = (
                            d,
                            ref_atom.name,
                            candidate_atom.name
                        )


        if best:

            d = best[0]

            if d <= STRONG_DISTANCE:
                score = 2
            else:
                score = 1


            if candidate_r not in active_network:

                active_network[candidate_r] = []


            active_network[candidate_r].append(
                (
                    d,
                    reference_r,
                    best[1],
                    best[2],
                    score
                )
            )


# ============================================================
# ANALYZE FE-BINDING RESIDUE CONNECTIONS
# ============================================================

fe_network = {}


for reference_r in FE_BINDING:

    reference_atoms = sidechain_atoms(
        reference_r
    )

    if not reference_atoms:
        continue


    for candidate_r in all_residues:

        # Don't analyze the 9 active-site residues themselves
        if candidate_r in ACTIVE_SITE:
            continue

        candidate_atoms = sidechain_atoms(
            candidate_r
        )

        if not candidate_atoms:
            continue


        best = None


        for ref_atom in reference_atoms:

            for candidate_atom in candidate_atoms:

                d = distance(
                    ref_atom,
                    candidate_atom
                )

                if d <= MAX_DISTANCE:

                    if (
                        best is None
                        or d < best[0]
                    ):

                        best = (
                            d,
                            reference_r,
                            ref_atom.name,
                            candidate_atom.name
                        )


        if best:

            d = best[0]

            if d <= STRONG_DISTANCE:
                score = 2
            else:
                score = 1


            if candidate_r not in fe_network:

                fe_network[candidate_r] = []


            fe_network[candidate_r].append(
                (
                    d,
                    reference_r,
                    best[2],
                    best[3],
                    score
                )
            )


# ============================================================
# ACTIVE-SITE NETWORK RESULTS
# ============================================================

print("\n")
print("=" * 80)
print("ACTIVE-SITE NETWORK")
print("=" * 80)

print(
    "Residue   Closest active-site residue   Distance   Score"
)

print("-" * 80)


for r, contacts in sorted(
    active_network.items(),
    key=lambda x: min(c[0] for c in x[1])
):

    contacts.sort(
        key=lambda x: x[0]
    )

    best = contacts[0]

    print(
        f"{residue_name(r)}{r:<7} "
        f"{residue_name(best[1])}{best[1]:<10} "
        f"{best[0]:>6.3f} A   "
        f"{best[4]}"
    )


# ============================================================
# Fe-BINDING NETWORK RESULTS
# ============================================================

print("\n")
print("=" * 80)
print("FE-BINDING RESIDUE NETWORK")
print("=" * 80)

print(
    "Residue   Closest Fe-binding residue   Distance   Score"
)

print("-" * 80)


for r, contacts in sorted(
    fe_network.items(),
    key=lambda x: min(c[0] for c in x[1])
):

    contacts.sort(
        key=lambda x: x[0]
    )

    best = contacts[0]

    print(
        f"{residue_name(r)}{r:<7} "
        f"{residue_name(best[1])}{best[1]:<10} "
        f"{best[0]:>6.3f} A   "
        f"{best[4]}"
    )


# ============================================================
# COMBINED NETWORK
# ============================================================

print("\n")
print("=" * 80)
print("COMBINED STRUCTURAL NETWORK")
print("=" * 80)

print(
    "Residue   Active-site evidence   Fe-binding evidence   Final"
)

print("-" * 80)


combined = {}


# Active-site evidence
for r, contacts in active_network.items():

    best = min(
        contacts,
        key=lambda x: x[0]
    )

    combined[r] = {
        "active_score": best[4],
        "active_distance": best[0],
        "active_reference": best[1],
        "fe_score": 0,
        "fe_distance": None,
        "fe_reference": None
    }


# Fe-binding evidence
for r, contacts in fe_network.items():

    best = min(
        contacts,
        key=lambda x: x[0]
    )

    if r not in combined:

        combined[r] = {
            "active_score": 0,
            "active_distance": None,
            "active_reference": None,
            "fe_score": best[4],
            "fe_distance": best[0],
            "fe_reference": best[1]
        }

    else:

        combined[r]["fe_score"] = best[4]
        combined[r]["fe_distance"] = best[0]
        combined[r]["fe_reference"] = best[1]


# ============================================================
# PRINT COMBINED RESULTS
# ============================================================

for r in sorted(combined):

    data = combined[r]

    final_score = max(
        data["active_score"],
        data["fe_score"]
    )


    if data["active_score"] > 0:

        active_text = (
            f"{data['active_score']} "
            f"({data['active_distance']:.2f} A "
            f"to {residue_name(data['active_reference'])}"
            f"{data['active_reference']})"
        )

    else:

        active_text = "0"


    if data["fe_score"] > 0:

        fe_text = (
            f"{data['fe_score']} "
            f"({data['fe_distance']:.2f} A "
            f"to {residue_name(data['fe_reference'])}"
            f"{data['fe_reference']})"
        )

    else:

        fe_text = "0"


    print(
        f"{residue_name(r)}{r:<7} "
        f"{active_text:<25} "
        f"{fe_text:<25} "
        f"{final_score}"
    )


# ============================================================
# SCORE 2
# ============================================================

SCORE_2 = []

for r, data in combined.items():

    if max(
        data["active_score"],
        data["fe_score"]
    ) == 2:

        SCORE_2.append(r)


print("\n")
print("=" * 80)
print("STRUCTURAL SCORE = 2")
print("=" * 80)

print(
    " ".join(
        f"{residue_name(r)}{r}"
        for r in sorted(SCORE_2)
    )
)


# ============================================================
# SCORE 1
# ============================================================

SCORE_1 = []

for r, data in combined.items():

    if max(
        data["active_score"],
        data["fe_score"]
    ) == 1:

        SCORE_1.append(r)


print("\n")
print("=" * 80)
print("STRUCTURAL SCORE = 1")
print("=" * 80)

print(
    " ".join(
        f"{residue_name(r)}{r}"
        for r in sorted(SCORE_1)
    )
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(
    "Maximum distance:",
    MAX_DISTANCE,
    "A"
)

print(
    "Number of active-site reference residues:",
    len(ACTIVE_SITE)
)

print(
    "Number of Fe-binding reference residues:",
    len(FE_BINDING)
)

print(
    "Structural score 2:",
    len(SCORE_2)
)

print(
    "Structural score 1:",
    len(SCORE_1)
)

print("=" * 80)

