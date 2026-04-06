"Generic head/tail splitter"

from pickle import NONE
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import pandas as pd
import os, statistics, argparse

# ---------------- CONFIG (defaults; can be overridden by CLI) ---------------- #

INPUT_EXCEL = None
SHEET = None
SMILES_COLUMN = None
DEFAULT_SMILES = ["CCCCCCCCCCCC(=O)NCCC[N+](C)(C)C(O)C[S](O)(=O)=O"]

MIN_TAIL_CARBONS = 3          # minimum carbon cluster size to qualify as tail
MAX_POLAR_DISTANCE = 1        # BFS layers to expand head from hetero/charged seeds
REFINE_WITH_CHARGES = True
CHARGE_PROMOTION_FACTOR = 2.0  # factor over median |q| of tail
MAX_PROMOTION_DISTANCE = 1

DEMOTE_EARLY_POLAR_SUBCLUSTER = True  # Demote early weak polar sub-cluster separated by neutral bridge from dominant polar cluster
BRIDGE_CHARGE_THRESH = 0.10  # |q| threshold to consider a bridge carbon weakly polar

DOMINANT_CLUSTER_CHARGE_RATIO = 2.0  # Minimum dominant/subcluster |q| sum ratio to allow demotion
REASSIGN_SMALL_TAILS = False     # Reassign very small tails into head
SMALL_TAIL_MAX_SIZE = 3          # Max size of a tail sub-component; if <=, move to head

# Global set of halogens (used to avoid classifying them as head and move them to tail)
HALOGENS = {9,17,35,53}  # F, Cl, Br, I

# Raw material inference thresholds
MIN_ETHO_UNITS_HEAD = 3        # minimum units to consider ethoxylated
MIN_SUGAR_OH = 4               # hetero (O) in head for sugar/polyol
MIN_POLYOL_OH = 3

# ------------- Basic utilities ------------- #
def mol_ok(smi):
    return Chem.MolFromSmiles(smi)

def is_hetero(atom):
    """Define hetero atoms for head seeds.
    We exclude halogens (F, Cl, Br, I) to avoid splitting perfluorinated chains; they are considered hydrophobic.
    """
    return atom.GetAtomicNum() not in (1,6,9,17,35,53)

def carbon_degree_carbon(atom):
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum()==6)

def hetero_neighbors(atom):
    return sum(1 for n in atom.GetNeighbors() if is_hetero(n))

def compute_gasteiger_abs(mol):
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception:
        return {i:0.0 for i in range(mol.GetNumAtoms())}
    charges={}
    for i in range(mol.GetNumAtoms()):
        a=mol.GetAtomWithIdx(i)
        try:
            q=a.GetDoubleProp('_GasteigerCharge')
        except Exception:
            q=0.0
        charges[i]=abs(q) if q is not None else 0.0
    return charges

def charge_refinement_details(mol, head, tails):
    """Compute the same promotion logic as refine_with_charges, but keep per-atom diagnostics."""
    details = {
        'enabled': REFINE_WITH_CHARGES and bool(tails),
        'charges': {},
        'near': set(),
        'component_thresholds': [],
        'promoted': [],
    }
    if not REFINE_WITH_CHARGES or not tails:
        return details

    charges = compute_gasteiger_abs(mol)
    details['charges'] = charges
    if not any(charges.values()):
        return details

    head_set = set(head)
    frontier = list(head_set)
    dist = {i: 0 for i in head_set}
    near = set(head_set)
    while frontier:
        idx = frontier.pop(0)
        if dist[idx] >= MAX_PROMOTION_DISTANCE:
            continue
        a = mol.GetAtomWithIdx(idx)
        for nb in a.GetNeighbors():
            nidx = nb.GetIdx()
            if nidx not in dist:
                dist[nidx] = dist[idx] + 1
                if dist[nidx] <= MAX_PROMOTION_DISTANCE:
                    near.add(nidx)
                    frontier.append(nidx)
    details['near'] = near

    for comp_id, comp in enumerate(tails, start=1):
        comp = set(comp)
        comp_ch = [charges[i] for i in comp]
        med = statistics.median(comp_ch) if comp_ch else 0.0
        thresh = med * CHARGE_PROMOTION_FACTOR if med > 0 else max(comp_ch or [0]) * 0.6
        details['component_thresholds'].append({
            'component_id': comp_id,
            'size': len(comp),
            'median_abs_charge': med,
            'threshold': thresh,
        })
        for idx in sorted(comp):
            atom = mol.GetAtomWithIdx(idx)
            hetero_nb = hetero_neighbors(atom)
            qualifies = idx in near and charges[idx] >= thresh and hetero_nb > 0
            if qualifies:
                details['promoted'].append({
                    'atom_idx': idx,
                    'symbol': atom.GetSymbol(),
                    'abs_charge': charges[idx],
                    'threshold': thresh,
                    'hetero_neighbors': hetero_nb,
                    'component_id': comp_id,
                    'distance_to_head': dist.get(idx),
                })
    return details

# ------------- Tail identification ------------- #
def hydrophobic_carbon_candidates(mol):
    cand=[]
    for a in mol.GetAtoms():
        if a.GetAtomicNum()!=6: continue
        # if hetero_neighbors(a) <= MAX_HETERO_NEIGHBORS_TAIL:
        cand.append(a.GetIdx())
    return cand

def carbon_components_from_candidates(mol, candidates):
    cand_set=set(candidates)
    comps=[]; seen=set()
    for idx in candidates:
        if idx in seen: continue
        stack=[idx]; comp=set()
        while stack:
            i=stack.pop()
            if i in comp: continue
            if i not in cand_set: continue
            comp.add(i)
            a=mol.GetAtomWithIdx(i)
            for nb in a.GetNeighbors():
                if nb.GetAtomicNum()==6 and nb.GetIdx() in cand_set and nb.GetIdx() not in comp:
                    stack.append(nb.GetIdx())
        seen|=comp
        comps.append(comp)
    return comps

def select_tail_components(mol, comps):
    if not comps:
        return []
    sizes=[len(c) for c in comps]
    largest=max(sizes)
    selected=[c for c in comps if len(c)>=MIN_TAIL_CARBONS and len(c)>=0.5*largest]
    if not selected:
        # fallback: take the largest if it exceeds minimum-2 (slightly relaxed)
        biggest=max(comps, key=len)
        if len(biggest)>=max(4, MIN_TAIL_CARBONS-2):
            selected=[biggest]
    return selected #[:MAX_TAILS_OUTPUT]

    # ------------- Head identification ------------- #
def head_seed_atoms(mol):
    seeds=set()
    for a in mol.GetAtoms():
        if is_hetero(a) or a.GetFormalCharge()!=0:
            seeds.add(a.GetIdx())
        else:
            if hetero_neighbors(a)>=2:
                seeds.add(a.GetIdx())
    return seeds

def expand_head(mol, seeds, tails):
    tail_atoms = set().union(*tails) if tails else set()
    head=set(seeds)
    queue=list(seeds)
    dist={i:0 for i in seeds}
    while queue:
        idx=queue.pop(0)
        a=mol.GetAtomWithIdx(idx)
        for nb in a.GetNeighbors():
            nidx=nb.GetIdx()
            if nidx in head or nidx in tail_atoms: continue
            # expand if hetero or within MAX_POLAR_DISTANCE proximity
            # Allow isolated halogens (not used as seeds) near the polar region; per-halogenated chains are corrected later.
            if is_hetero(nb) or dist[idx] < MAX_POLAR_DISTANCE:
                head.add(nidx)
                dist[nidx]=dist[idx]+1
                queue.append(nidx)
    if not head:  # fallback: at least one hetero atom if present
        for a in mol.GetAtoms():
            if is_hetero(a):
                head.add(a.GetIdx())
                break
    return head

# ------------- Charge-based refinement ------------- #
def refine_with_charges(mol, head, tails):
    if not REFINE_WITH_CHARGES or not tails:
        return head, tails, False
    details = charge_refinement_details(mol, head, tails)
    charges = details['charges']
    if not any(charges.values()):
        return head, tails, False
    head_set=set(head)
    promoted=set()
    near = details['near']
    threshold_map = {d['component_id']: d['threshold'] for d in details['component_thresholds']}
    refined=[]
    for comp_id, comp in enumerate(tails, start=1):
        if not comp:
            refined.append(comp); continue
        thresh = threshold_map.get(comp_id, 0.0)
        new_comp=set(comp)
        for idx in list(comp):
            if idx in near and charges[idx] >= thresh and hetero_neighbors(mol.GetAtomWithIdx(idx))>0:
                promoted.add(idx)
                new_comp.discard(idx)
        refined.append(new_comp)
    if promoted:
        head_set |= promoted
    return head_set, refined, bool(promoted)

def diagnose_charge_refinement(smiles_list):
    rows = []
    for smi in smiles_list:
        mol = mol_ok(smi)
        if not mol:
            rows.append({'Surfactant': smi, 'Valid': False})
            continue
        candidates = hydrophobic_carbon_candidates(mol)
        comps = carbon_components_from_candidates(mol, candidates)
        tails = select_tail_components(mol, comps)
        seeds = head_seed_atoms(mol)
        head = expand_head(mol, seeds, tails)
        cleaned = [set(a for a in t if a not in head) for t in tails]
        details = charge_refinement_details(mol, head, cleaned)

        if not details['promoted']:
            rows.append({
                'Surfactant': smi,
                'Valid': True,
                'ChargeRefined': False,
                'PromotedAtomIdx': None,
                'AtomSymbol': None,
                'AbsCharge': None,
                'Threshold': None,
                'HeteroNeighbors': None,
                'ComponentId': None,
                'Reason': 'No atom met near-head + charge-threshold + hetero-neighbor criteria'
            })
            continue

        for item in details['promoted']:
            rows.append({
                'Surfactant': smi,
                'Valid': True,
                'ChargeRefined': True,
                'PromotedAtomIdx': item['atom_idx'],
                'AtomSymbol': item['symbol'],
                'AbsCharge': item['abs_charge'],
                'Threshold': item['threshold'],
                'HeteroNeighbors': item['hetero_neighbors'],
                'ComponentId': item['component_id'],
                'Reason': f"near head, |q| >= threshold, hetero_neighbors={item['hetero_neighbors']}"
            })
    return pd.DataFrame(rows)

# ------------- Assembly and HI calculation ------------- #
def frag_smiles(mol, atoms):
    if not atoms: return ''
    return Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(atoms), canonical=True)  # type: ignore[attr-defined]

def compute_hi(mol, head_smiles):
    hm = Chem.MolFromSmiles(head_smiles)
    if not hm:
        return None
    return 20 * (Descriptors.MolWt(hm)/Descriptors.MolWt(mol))

# ------------- Main per-molecule pipeline ------------- #
EXCLUDE_COUNTERIONS = True  # Exclude small disconnected charged fragments (counterions)

def split_head_tail(mol):
    # 1. Preliminary tails
    candidates = hydrophobic_carbon_candidates(mol)
    comps = carbon_components_from_candidates(mol, candidates)
    tails = select_tail_components(mol, comps)
    # 2. Head
    seeds = head_seed_atoms(mol)
    head = expand_head(mol, seeds, tails) 
    counterion_atoms = set()  
    removed_counterion_atoms = set()
    # 2a. Counterions elimination
    counterion_smiles_list = []
    counterion_removed_flag = False
    if EXCLUDE_COUNTERIONS:
        try:
            frags = Chem.GetMolFrags(mol)  
        except Exception:
            frags = []
        if frags and len(frags) > 1:
            def heavy_count(f):
                return sum(1 for idx in f if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
            main_frag = max(frags, key=lambda f: (heavy_count(f), len(f)))
            main_set = set(main_frag)
            for frag in frags:
                if frag is main_frag:
                    continue
                fset = set(frag)
                heavy = [mol.GetAtomWithIdx(i) for i in fset if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
                charge_total = sum(a.GetFormalCharge() for a in heavy)
                # criteria: total charge non 0 and small size (<= 5 heavy atoms)
                if charge_total != 0 and len(heavy) <= 5:
                    counterion_atoms |= fset
                    frag_smiles2 = Chem.MolFragmentToSmiles(mol, list(fset), isomericSmiles=True)
                    counterion_smiles_list.append(frag_smiles2)
            if counterion_atoms:
                head -= counterion_atoms
                for t in tails:
                    t -= counterion_atoms
                removed_counterion_atoms = set(counterion_atoms)
                counterion_removed_flag = True
    # 3. Clean overlap (if any tail shares head atoms, remove those atoms from the tail)
    cleaned=[]
    for t in tails:
        cleaned.append(set(a for a in t if a not in head))
    tails = cleaned
    # 3a. Promote rings with hetero or charged atoms entirely to head (e.g., pyridinium)
    head, tails = promote_hetero_charged_rings(mol, head, tails)
    # 3a2. Promote carbohydrate rings (glucose-like) entirely to head
    head, tails = promote_carbohydrate_rings(mol, head, tails)
    # 3a3. Promote open poly-hydroxylated chains (linear polyols)
    head, tails = promote_polyol_chains(mol, head, tails)
    
    # 3b. Force full benzene rings into tail (avoid splitting the same phenyl between head/tail)
    if tails:
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if len(ring)==6 and all(mol.GetAtomWithIdx(i).GetAtomicNum()==6 and mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                rset = set(ring)
                for t in tails:
                    if t & rset:  
                        if not rset.issubset(t):
                            t |= rset
                            head -= rset
                        break 
    # 4. Charge-based refinement
    head, tails, charge_ref = refine_with_charges(mol, head, tails)
    # 5. Boundary adjustments (later steps may adjust linkers)
    boundary_adjusted = False
    sugar_like_head = sum(1 for i in head if mol.GetAtomWithIdx(i).GetAtomicNum()==8) >= MIN_SUGAR_OH
    # 5b. Enforce full benzene integrity after adjustments
    head, tails = enforce_full_benzenes(mol, head, tails)
    # 5c. Re-promote linking carbonyl if sugar-like head
    if sugar_like_head:
        head, tails = promote_linking_amide_carbonyl(mol, head, tails)
    # 5d. Additional rule: no purely carbon aromatic ring (no hetero/charge) should remain isolated in head.
    #     If a ring is in head (partially or fully), force entire ring into a tail (existing or new).
    if True:
        ri = mol.GetRingInfo()
        ring_atoms_list = ri.AtomRings()
        for ring in ring_atoms_list:
            # consider only pure C aromatic rings
            if not all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum()==6 for i in ring):
                continue
            ring_set = set(ring)
            # exclude if already promoted as hetero/charged (not applicable here because these are pure C)
            # or if the ring is already fully in a tail
            in_head = ring_set & head
            if not in_head:
                continue  # nothing in head -> already fully in tail or outside
            # move full ring to a tail
            target_tail = None
            for t in tails:
                if t & ring_set:
                    target_tail = t
                    break
            if target_tail is None:
                # choose largest tail if present, otherwise create a new one
                if tails:
                    target_tail = max(tails, key=len)
                else:
                    target_tail = set()
                    tails.append(target_tail)
            target_tail |= ring_set
            head -= ring_set
    # 6. Head refinement into multiple functional fragments (e.g., multiple sulfates)
    head, head_fragments, tails = refine_heads_split(mol, head, tails)
    # 6b2. Additional failsafe: if head becomes empty, rescue a hetero or first carbon
    head_failsafe=False
    if not head:
        hetero=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1,6)]
        if hetero:
            h=hetero[0]
            head.add(h)
            for t in tails:
                if h in t: t.remove(h)
            # Add a neighboring carbon if exists
            for nb in mol.GetAtomWithIdx(h).GetNeighbors():
                if nb.GetAtomicNum()==6:
                    head.add(nb.GetIdx())
                    for t in tails:
                        if nb.GetIdx() in t: t.remove(nb.GetIdx())
                    break
        else:
            # Take atom 0 as minimal head
            if mol.GetNumAtoms()>0:
                head.add(0)
                for t in tails:
                    if 0 in t: t.remove(0)
        head_failsafe=True
    # 6b3. Optional expansion: if only O as hetero (simple alcohols/polyols), promote all O atoms and their alpha carbons
    def expand_simple_alcohol_head(mol, head_set, tails_list):
        atoms = list(mol.GetAtoms())
        # Detect if there are hetero atoms other than O
        has_other_hetero = any(a.GetAtomicNum() not in (1,6,8) for a in atoms)
        if has_other_hetero:
            return head_set, False
        o_indices=[a.GetIdx() for a in atoms if a.GetAtomicNum()==8]
        if not o_indices:
            return head_set, False
        changed=False
        for oidx in o_indices:
            head_set.add(oidx)
            oatom=mol.GetAtomWithIdx(oidx)
            for nb in oatom.GetNeighbors():
                if nb.GetAtomicNum()==6:
                    cid=nb.GetIdx()
                    if cid not in head_set:
                        head_set.add(cid); changed=True
                        # remove from tails if present
                    for t in tails_list:
                        if cid in t: t.remove(cid)
                    # remove oxygen from tails
            for t in tails_list:
                if oidx in t: t.remove(oidx)
        return head_set, changed
    head, expanded_alcohol = expand_simple_alcohol_head(mol, head, tails)
    if expanded_alcohol:
        head_failsafe = True
    # 6c. Failsafe: no unassigned atoms
    all_atoms = set(range(mol.GetNumAtoms()))
    tail_union = set().union(*tails) if tails else set()
    unassigned = all_atoms - head - tail_union
    if EXCLUDE_COUNTERIONS and removed_counterion_atoms:
        unassigned -= removed_counterion_atoms
    if unassigned:
        # Rule: simply add them to head to avoid losing atoms
        head |= unassigned
    # 6d. Specific cleanup: remove aromatic "stub" (isolated aromatic carbon) remaining in head adjacent to O.
    #     Observed case: Head = cOCCO. Rule: any aromatic carbon in head whose full ring is not in head,
    #     and most of the ring is in tail, gets moved to tail.
    if head:
        ri = mol.GetRingInfo()
        ring_atoms_list = ri.AtomRings()
        for ring in ring_atoms_list:
            ring_set = set(ring)
            # only consider aromatic rings made entirely of C
            if not all(mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum()==6 for i in ring_set):
                continue
            head_in_ring = ring_set & head
            if 0 < len(head_in_ring) < len(ring_set):
                # fragmented; evaluate if it is a stub (1 or 2 carbons) attached to external hetero
                if len(head_in_ring) <= 2:
                    # move those carbons to a tail (prefer a tail that already contains part of the ring, or the largest)
                    target_tail = None
                    for t in tails:
                        if t & ring_set:
                            target_tail = t; break
                    if target_tail is None:
                        if tails:
                            target_tail = max(tails, key=len)
                        else:
                            target_tail = set(); tails.append(target_tail)
                    target_tail |= head_in_ring
                    head -= head_in_ring
        # Second pass: if head becomes empty due to this cleanup, rescue the nearest oxygen (if any)
        if not head:
            oxy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==8]
            if oxy:
                oidx = oxy[0]; head.add(oidx)
                # Remove from tail
                for t in tails:
                    if oidx in t: t.remove(oidx)
    # 6e. Promote fully per-halogenated chains (perfluoro, perchloro, etc.) to tail if partially in head.
    # Treat C chains enriched in halogens (F, Cl, Br, I) as hydrophobic if heavily halogenated.
    # Move the carbon component and ALL directly bound halogens.
    # Criteria:
    #   - Carbon component where non-carbon neighbors are only halogens or (on at most one carbon) a terminal O/S.
    #   - Fraction of carbons with ≥1 halogen ≥ 70%.Modified to 100%
    #   - Intersects head (part was misclassified as head).
    #   - Does not completely remove real hetero head (non-halogen) leaving it empty.
    perfluoro_moved = False  # mantiene nombre legado
    perhalogen_moved = False
    visited_ph=set()
    for a in mol.GetAtoms():
        idx=a.GetIdx()
        if idx in visited_ph: continue
        if a.GetAtomicNum()!=6: continue
        # BFS of candidate carbon component
        stack=[idx]; comp=set(); qualifies=True
        carbon_with_SO_neighbors=set()
        while stack:
            i=stack.pop()
            if i in comp: continue
            ai=mol.GetAtomWithIdx(i)
            if ai.GetAtomicNum()!=6:
                qualifies=False; continue
            comp.add(i); visited_ph.add(i)
            for nb in ai.GetNeighbors():
                nZ=nb.GetAtomicNum()
                if nZ not in (1,6,*HALOGENS,8,16):  # only C, H, halogens and terminal O/S
                    qualifies=False
                if nZ in (8,16):
                    carbon_with_SO_neighbors.add(i)
                if nZ==6 and nb.GetIdx() not in comp:
                    stack.append(nb.GetIdx())
        if not qualifies or not comp:
            continue
        if len(carbon_with_SO_neighbors) > 1:
            continue
        # Fraction of halogenated carbons
        halogenated_carbons=0
        for ci in comp:
            if any(nb.GetAtomicNum() in HALOGENS for nb in mol.GetAtomWithIdx(ci).GetNeighbors()):
                halogenated_carbons+=1
        if halogenated_carbons==0 or halogenated_carbons/len(comp) < 1:
            continue
        # Does it intersect head?
        if not (comp & head):
            continue
        head_remaining = head - comp
        has_real_hetero = any(mol.GetAtomWithIdx(h).GetAtomicNum() not in (1,6,*HALOGENS) for h in head)
        if not head_remaining and not has_real_hetero:
            # do not move: keep a minimal polar head
            continue
        # Gather halogens bound to the component
        halogen_atoms=set()
        for ci in comp:
            for nb in mol.GetAtomWithIdx(ci).GetNeighbors():
                if nb.GetAtomicNum() in HALOGENS:
                    halogen_atoms.add(nb.GetIdx())
        # Move to tail
        target_tail=None
        for t in tails:
            if t & comp:
                target_tail=t; break
        if target_tail is None:
            if tails:
                target_tail=max(tails, key=len)
            else:
                target_tail=set(); tails.append(target_tail)
        target_tail |= comp | halogen_atoms
        head -= (comp | halogen_atoms)
        perhalogen_moved = True
        if any(mol.GetAtomWithIdx(h).GetAtomicNum()==9 for h in halogen_atoms):
            perfluoro_moved = True
    # 6e2. Generic post-pass: any halogen remaining in head whose non-H neighbors are all in tail is moved to tail.
    halogen_post_moved = False
    if tails:
        tail_union_all = set().union(*tails)
        for hidx in list(head):
            a = mol.GetAtomWithIdx(hidx)
            if a.GetAtomicNum() not in HALOGENS:
                continue
            # Ignore if it has any real hetero neighbor (O,N,S, etc.) in head (possible part of polar functional group)
            neighs = [nb.GetIdx() for nb in a.GetNeighbors() if nb.GetAtomicNum()!=1]
            if not neighs:
                continue
            # If all non-H neighbors are in tail (and none in head except the halogen itself) => move
            if all((nid in tail_union_all) for nid in neighs):
                # Place in the tail that contains the first neighbor
                placed=False
                for t in tails:
                    if neighs[0] in t:
                        t.add(hidx); placed=True; break
                if not placed:
                    # Fallback: add to the largest tail
                    max(tails, key=len).add(hidx)
                head.discard(hidx)
                halogen_post_moved = True
                if a.GetAtomicNum()==9:
                    perfluoro_moved = True

    carboxylate_repromoted = False
    for t in tails:
        for cidx in list(t):
            aC = mol.GetAtomWithIdx(cidx)
            if aC.GetAtomicNum()!=6:
                continue
            dbl_O = None
            single_O = None
            other_carb_tail = False
            for b in aC.GetBonds():
                other = b.GetOtherAtom(aC)
                if b.GetBondType().name == 'DOUBLE' and other.GetAtomicNum()==8:
                    dbl_O = other.GetIdx()
                elif b.GetBondType().name == 'SINGLE' and other.GetAtomicNum()==8:
                    single_O = other.GetIdx()
                elif other.GetAtomicNum()==6 and other.GetIdx() in t:
                    other_carb_tail = True
            if not (dbl_O and single_O and other_carb_tail):
                continue
            o1 = mol.GetAtomWithIdx(dbl_O)
            o2 = mol.GetAtomWithIdx(single_O)

            polar_flag = (o2.GetFormalCharge() < 0) or (o1.GetFormalCharge() != 0) or (o2.GetFormalCharge() != 0)
            if polar_flag or single_O not in t:
                moved=False
                for atom_id in (cidx, dbl_O, single_O):
                    if atom_id in t and atom_id not in head:
                        t.discard(atom_id)
                        head.add(atom_id)
                        moved=True
                if moved:
                    carboxylate_repromoted = True

    
    # 6g. Ethoxylate boundary trim: exclude the first alkyl carbon immediately after the last EO oxygen (O-CH2-CH2-R -> move R-CH2 alpha)
    ethox_boundary_trimmed = False
    try:
        # Heuristic: find head oxygens (O in head) with exactly 2 carbon neighbors both in head.
        # Identify terminal EO oxygen: one side leads (in ≤2 steps) to another oxygen within head (internal EO chain) and the other side leads
        # to a hydrocarbon region without nearby oxygens (≥2 steps without finding O). That hydrocarbon carbon is reclassified to tail.
        head_oxygens = [i for i in head if mol.GetAtomWithIdx(i).GetAtomicNum()==8]
        oxy_to_move_carbons = []
        for oidx in head_oxygens:
            oatom = mol.GetAtomWithIdx(oidx)
            neigh_carbons = [nb.GetIdx() for nb in oatom.GetNeighbors() if nb.GetAtomicNum()==6 and nb.GetIdx() in head]
            if len(neigh_carbons)!=2:
                continue
            c1, c2 = neigh_carbons
            # Detect another oxygen within ≤2 steps moving away from oidx through a specific carbon
            def oxygen_ahead(start_c):
                visited={oidx}
                frontier=[(start_c,0)]
                while frontier:
                    cur,depth=frontier.pop(0)
                    if depth>2: break
                    if cur!=start_c and mol.GetAtomWithIdx(cur).GetAtomicNum()==8 and cur in head:
                        return True
                    a_cur = mol.GetAtomWithIdx(cur)
                    for nb in a_cur.GetNeighbors():
                        nid=nb.GetIdx()
                        if nid in visited: continue
                        # do not return to the initial oxygen except for the implicit first step
                        visited.add(nid)
                        # allow traversal through carbons or oxygens in head
                        if nb.GetAtomicNum() in (6,8):
                            frontier.append((nid, depth+1))
                return False
            c1_internal = oxygen_ahead(c1)
            c2_internal = oxygen_ahead(c2)
            # Exactly one must be internal EO; the other is tail candidate
            if c1_internal == c2_internal:
                continue
            tail_candidate = c1 if not c1_internal else c2
            # Confirm tail_candidate has no other oxygen neighbors (only this oidx) and is aliphatic (no extra hetero)
            tc_atom = mol.GetAtomWithIdx(tail_candidate)
            hetero_non_oxygen = any(nb.GetAtomicNum() not in (1,6,8) for nb in tc_atom.GetNeighbors())
            oxygen_neighbors = [nb for nb in tc_atom.GetNeighbors() if nb.GetAtomicNum()==8]
            if hetero_non_oxygen or len(oxygen_neighbors)!=1:
                continue
            # Check that beyond (2-3 steps) there is no oxygen -> hydrophobic chain
            def oxygen_within(start, max_depth=3):
                visited={oidx}
                frontier=[(start,0)]
                while frontier:
                    cur,depth=frontier.pop(0)
                    if depth>max_depth: continue
                    if cur!=tail_candidate and mol.GetAtomWithIdx(cur).GetAtomicNum()==8 and cur in head:
                        return True
                    a_cur=mol.GetAtomWithIdx(cur)
                    for nb in a_cur.GetNeighbors():
                        nid=nb.GetIdx()
                        if nid in visited: continue
                        visited.add(nid)
                        if mol.GetAtomWithIdx(nid).GetAtomicNum() in (6,8):
                            frontier.append((nid, depth+1))
                return False
            if oxygen_within(tail_candidate):
                # Found another nearby oxygen -> likely still part of internal EO
                continue
            oxy_to_move_carbons.append(tail_candidate)
        if oxy_to_move_carbons:
            for cidx in oxy_to_move_carbons:
                if cidx not in head:
                    continue
                # Choose destination tail: one that already connects to any neighboring carbons; if none, largest existing tail
                neighbors = {nb.GetIdx() for nb in mol.GetAtomWithIdx(cidx).GetNeighbors() if nb.GetAtomicNum()==6}
                placed=False
                for t in tails:
                    if t & neighbors:
                        t.add(cidx); placed=True; break
                if not placed:
                    if tails:
                        max(tails, key=len).add(cidx)
                    else:
                        # Create new tail
                        tails.append({cidx})
                head.discard(cidx)
                ethox_boundary_trimmed = True
    except Exception:
        pass

    # 6l-4. Demotion of early weak polar sub-cluster separated by neutral bridge from dominant cluster
    early_polar_cluster_demoted = False
    if DEMOTE_EARLY_POLAR_SUBCLUSTER and head and len(head_fragments) > 1:
        gasteiger_charges = compute_gasteiger_abs(mol)
        comp_metrics = []  # (fragment_set, abs_charge_sum, formal_charge_sum, hetero_count)
        for frag in head_fragments:
            abs_sum = sum(gasteiger_charges.get(i,0.0) for i in frag)
            formal_sum = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in frag)
            hetero_count = sum(1 for i in frag if (mol.GetAtomWithIdx(i).GetAtomicNum() not in (1,6,*HALOGENS)))
            comp_metrics.append((frag, abs_sum, formal_sum, hetero_count))
        def dominance_key(item):
            frag, abs_sum, formal_sum, hetero_count = item
            return (
                1 if formal_sum!=0 else 0,
                hetero_count,
                abs_sum
            )
        dominant = max(comp_metrics, key=dominance_key)
        dom_frag, dom_abs, dom_formal, dom_hetero = dominant
        # Candidates: no formal charge and significantly lower abs_sum
        candidates = [m for m in comp_metrics if m[0] is not dom_frag]

        # Prebuild quick head access for BFS
        head_set_local = set(head)

        from collections import deque
        def shortest_path_between_sets(starts, targets, allowed):
            target_set = set(targets)
            visited = set(starts)
            parent = {}
            dq = deque(starts)
            while dq:
                cur = dq.popleft()
                if cur in target_set:
                    # reconstruct
                    path = [cur]
                    while cur in parent:
                        cur = parent[cur]
                        path.append(cur)
                    path.reverse()
                    return path
                a_cur = mol.GetAtomWithIdx(cur)
                for nb in a_cur.GetNeighbors():
                    nid = nb.GetIdx()
                    if nid not in allowed or nid in visited:
                        continue
                    visited.add(nid)
                    parent[nid] = cur
                    dq.append(nid)
            return None

        demote_total = set()
        for frag, abs_sum, formal_sum, hetero_count in candidates:
            if formal_sum != 0:
                continue  # has formal charge, potential real head (avoid false gemini detection)
            if dom_abs <= 0 or abs_sum <= 0:
                continue
            ratio = dom_abs / max(abs_sum, 1e-6)
            if ratio < DOMINANT_CLUSTER_CHARGE_RATIO:
                continue  # dominant cluster is not sufficiently more polar
            # Find shortest path in the head subgraph
            path = shortest_path_between_sets(list(frag), list(dom_frag), head_set_local)
            if not path:
                continue
            # Exclude terminal nodes belonging to frag or dom_frag to evaluate bridge
            inner = [n for n in path if (n not in frag and n not in dom_frag)]
            if not inner:
                continue
            # Count consecutive weakly polar carbons
            seq_count = 0
            max_seq = 0
            for n in inner:
                a = mol.GetAtomWithIdx(n)
                if a.GetAtomicNum()==6 and gasteiger_charges.get(n,0.0) < BRIDGE_CHARGE_THRESH:
                    seq_count += 1
                    max_seq = max(max_seq, seq_count)
                else:
                    seq_count = 0
            demote_frag = set(frag)
            bridge_carbons = {n for n in inner if mol.GetAtomWithIdx(n).GetAtomicNum()==6 and gasteiger_charges.get(n,0.0) < BRIDGE_CHARGE_THRESH}
            demote_set = demote_frag | bridge_carbons
            demote_total |= demote_set

        if demote_total:
            # Select destination tail (adjacent or largest)
            for u in demote_total:
                target = None
                ua = mol.GetAtomWithIdx(u)
                nbs = {nb.GetIdx() for nb in ua.GetNeighbors()}
                for t in tails:
                    if t & nbs:
                        target = t; break
                if target is None:
                    target = max(tails, key=len) if tails else None
                if target is None:
                    target = set(); tails.append(target)
                if u in head:
                    head.remove(u)
                    target.add(u)
                    early_polar_cluster_demoted = True
            if early_polar_cluster_demoted:
                # Recompute final head fragments
                def _recompute_head_fragments2(hset):
                    comps = []
                    visited = set()
                    for idx in hset:
                        if idx in visited: continue
                        stack=[idx]; comp=set()
                        while stack:
                            x=stack.pop()
                            if x in visited or x not in hset: continue
                            visited.add(x); comp.add(x)
                            ax = mol.GetAtomWithIdx(x)
                            for nb in ax.GetNeighbors():
                                nid = nb.GetIdx()
                                if nid in hset and nid not in visited:
                                    stack.append(nid)
                        comps.append(comp)
                    return comps
                head_fragments = _recompute_head_fragments2(head)
    # 6l-5. Reassignment of small tails: no tail cluster of size <= SMALL_TAIL_MAX_SIZE should remain isolated as tail
    small_tail_components_reassigned = False
    if REASSIGN_SMALL_TAILS and tails:
        # Iterate each tail; if a full tail is small, move it.
        # Also decompose each tail into connected subcomponents for robustness.
        new_tails = []
        subs_small_moved = []  # store moved small subcomponents in case tail rescue is needed
        for t in tails:
            # Build connected subcomponents within t
            subs = []
            seen_sub = set()
            for aidx in list(t):
                if aidx in seen_sub: continue
                stack=[aidx]; comp=set()
                while stack:
                    x=stack.pop()
                    if x in comp: continue
                    if x not in t: continue
                    comp.add(x); seen_sub.add(x)
                    ax = mol.GetAtomWithIdx(x)
                    for nb in ax.GetNeighbors():
                        nid = nb.GetIdx()
                        if nid in t and nid not in comp:
                            stack.append(nid)
                subs.append(comp)
            for comp in subs:
                if len(comp) <= SMALL_TAIL_MAX_SIZE:
                    # Move to head (mark for possible rescue)
                    subs_small_moved.append(comp)
                    for ai in comp:
                        head.add(ai)
                        small_tail_components_reassigned = True
                else:
                    new_tails.append(comp)
        tails = [set(c) for c in new_tails]
        # If all tails were small and reassigned, rescue the largest to ensure at least one tail
        if not tails and subs_small_moved:
            largest_small = max(subs_small_moved, key=len)
            # remove from head and restore as tail
            for ai in largest_small:
                if ai in head:
                    head.remove(ai)
            tails = [set(largest_small)]
            # Keep small_tail_components_reassigned=True (small ones reassigned except the rescued)
        if small_tail_components_reassigned:
            # Recompute head fragments after incorporating/rescuing these atoms
            def _recompute_head_fragments_small(hset):
                comps=[]; visited=set()
                for idx in hset:
                    if idx in visited: continue
                    st=[idx]; comp=set()
                    while st:
                        j=st.pop()
                        if j in comp: continue
                        if j not in hset: continue
                        comp.add(j); visited.add(j)
                        aj=mol.GetAtomWithIdx(j)
                        for nb in aj.GetNeighbors():
                            nid=nb.GetIdx()
                            if nid in hset and nid not in comp:
                                st.append(nid)
                    comps.append(comp)
                return comps
            head_fragments = _recompute_head_fragments_small(head)

    if EXCLUDE_COUNTERIONS and removed_counterion_atoms:
        head -= removed_counterion_atoms  # final purge
    head_smiles = frag_smiles(mol, head)
    head_frag_smiles = [frag_smiles(mol, h) for h in head_fragments]
    tail_smiles_list = [frag_smiles(mol, t) for t in tails]
    combined = set().union(*tails) if tails else set()
    combined_smiles = frag_smiles(mol, combined)
    combined_list= [combined]
    hi = compute_hi(mol, head_smiles)
    return {
        'Head': head_smiles,
        'HeadAtoms': head,
        'HeadAtomsList': head_fragments,
        'HeadCount': len(head_fragments),
        'HeadList': head_frag_smiles[:3],
        'TailCombined': combined_smiles,
        'TailAtomsList': combined_list,
        'TailCount': len(tail_smiles_list),
        'TailList': tail_smiles_list,
        'HI': hi,
        'CounterionList': counterion_smiles_list,
        'BoundaryAdjusted': boundary_adjusted,
        'HeadFailsafe': head_failsafe
        ,'PerfluoroMoved': perfluoro_moved
        ,'PerhalogenMoved': perhalogen_moved
        ,'HalogenPostMoved': halogen_post_moved
        # Note: isolated halogens may remain in head; per-halogenated chains are moved.
        ,'CarboxylateRepromoted': carboxylate_repromoted
        ,'EarlyPolarClusterDemoted': early_polar_cluster_demoted
        ,'SmallTailComponentsReassigned': small_tail_components_reassigned
    }

# ------------- Boundary adjustment to include acyl carbonyl C=O in tail ------------- #
def adjust_acyl_boundary(mol, head_set, tails_list):
    head_set = set(head_set)
    adjusted = False
    # Pre-index tails for quick lookup
    for t in tails_list:
        # Find head carbonyls adjacent to tail carbon and head hetero
        candidate_carbons = [idx for idx in list(head_set) if mol.GetAtomWithIdx(idx).GetAtomicNum()==6]
        for c_idx in candidate_carbons:
            a = mol.GetAtomWithIdx(c_idx)
            # Must have a double bond to O (carbonyl) and a bridging O/N (ester/amide link)
            dbl_O = None            # double-bond oxygen
            linking_O = None        # single-bond O (or other hetero) linking toward head
            hetero_bridge = False
            tail_contact = False
            for b in a.GetBonds():
                bt = b.GetBondType().name
                other = b.GetOtherAtom(a)
                o_idx = other.GetIdx()
                if bt == 'DOUBLE' and other.GetAtomicNum()==8:
                    dbl_O = o_idx
                elif bt == 'SINGLE':
                    # Hetero neighbor in head (O or N) linking toward polar part
                    if other.GetAtomicNum() in (7,8) and o_idx in head_set:
                        hetero_bridge = True
                        if other.GetAtomicNum()==8:  # prefer moving ester bridging O
                            linking_O = o_idx
                    # Carbon neighbor in tail
                    if other.GetIdx() in t and other.GetAtomicNum()==6:
                        tail_contact = True
            if dbl_O is not None and hetero_bridge and tail_contact:
                # Move full carbonyl: C, =O and bridging O (if any) to tail
                moved_any=False
                if c_idx in head_set:
                    head_set.remove(c_idx)
                    t.add(c_idx)
                    moved_any=True
                if dbl_O in head_set:
                    head_set.remove(dbl_O)
                    t.add(dbl_O)
                    moved_any=True
                if linking_O is not None and linking_O in head_set:
                    head_set.remove(linking_O)
                    t.add(linking_O)
                    moved_any=True
                if moved_any:
                    adjusted = True
    return head_set, tails_list, adjusted

# ------------- Re-promotion of linking amide carbonyl (amide/ester linker) ------------- #
def promote_linking_amide_carbonyl(mol, head_set, tails_list):
        # """Reincorporate the carbonyl C=O (C and O atoms) of an amide/ester linker into the head when connecting
        # a sugar/polyol head (many oxygens in head) to a tail containing the nitrogen (or alkyl chain).

        # Criteria:
        #     - Pattern C (with =O) + neighboring N (amide-like) where C or its =O are in tail but
        #         there is at least one neighbor (O or C) belonging to the oxygen-rich head.
        #     - Only runs if head has >= MIN_SUGAR_OH oxygens (sugar_like_head).
        #     - Moves carbonyl C and its carbonyl O from tail to head; does not move the nitrogen.
        # """
    head_set = set(head_set)
    oxy_head = sum(1 for i in head_set if mol.GetAtomWithIdx(i).GetAtomicNum()==8)
    if oxy_head < MIN_SUGAR_OH:
        return head_set, tails_list
    # Traverse tail atoms looking for amide/ester carbonyls
    for t in tails_list:
        for idx in list(t):
            a = mol.GetAtomWithIdx(idx)
            if a.GetAtomicNum()!=6:
                continue
            dbl_O = None
            neigh_N = None
            for b in a.GetBonds():
                other = b.GetOtherAtom(a)
                bt = b.GetBondType().name
                if bt == 'DOUBLE' and other.GetAtomicNum()==8:
                    dbl_O = other.GetIdx()
                elif other.GetAtomicNum()==7:
                    neigh_N = other.GetIdx()
            if dbl_O is None or neigh_N is None:
                continue
            # Connected to head? (any O or C neighbor in head)
            linked_head = (
                any(nb.GetIdx() in head_set and nb.GetAtomicNum()==8 for nb in a.GetNeighbors()) or
                any(nb.GetIdx() in head_set and nb.GetAtomicNum()==6 for nb in a.GetNeighbors())
            )
            if linked_head:
                if idx in t:
                    t.discard(idx)
                    head_set.add(idx)
                if dbl_O in t:
                    t.discard(dbl_O)
                    head_set.add(dbl_O)
    return head_set, tails_list

# ------------- Promotion of hetero/charged rings to head ------------- #
def promote_hetero_charged_rings(mol, head_set, tails_list):
    """Any aromatic ring containing at least one hetero atom (not C/H) or an atom with formal charge is moved entirely to head.
    Removes its atoms from any tails that contain them."""
    ri = mol.GetRingInfo()
    ring_atoms_list = ri.AtomRings()
    head_set = set(head_set)
    for ring in ring_atoms_list:
        ring_set = set(ring)
        has_hetero_or_charge = False
        for idx in ring:
            a = mol.GetAtomWithIdx(idx)
            if a.GetFormalCharge() != 0 or a.GetAtomicNum() not in (1,6):
                has_hetero_or_charge = True
                break
        if has_hetero_or_charge:
            # remove from any tails
            for t in tails_list:
                if t & ring_set:
                    t -= ring_set
            head_set |= ring_set
    return head_set, tails_list

# ------------- Promote small carbon bridges between head fragments ------------- #
def promote_head_bridges(mol, head_set, tails_list, max_bridge_len=2):
    """If small fragments (only C/H) of length <= max_bridge_len connect two head regions,
    add them to head. Avoid losing intermediate carbons like in -CO-CH2-CH2-NH- chains when both ends are already head.
    """
    head_set = set(head_set)
    tail_union = set().union(*tails_list) if tails_list else set()
    all_atoms = set(range(mol.GetNumAtoms()))
    candidate = all_atoms - head_set - tail_union
    if not candidate:
        return head_set, tails_list
    # build candidate components made only of carbons
    visited=set()
    bridges=[]
    for idx in list(candidate):
        if idx in visited:
            continue
        a=mol.GetAtomWithIdx(idx)
        if a.GetAtomicNum()!=6:
            continue
        stack=[idx]; comp=set()
        only_carbon=True
        while stack:
            i=stack.pop()
            if i in comp: continue
            ai=mol.GetAtomWithIdx(i)
            if ai.GetAtomicNum()!=6: only_carbon=False
            comp.add(i)
            visited.add(i)
            for nb in ai.GetNeighbors():
                j=nb.GetIdx()
                if j in head_set or j in tail_union or j in comp:
                    continue
                stack.append(j)
        if only_carbon and 0 < len(comp) <= max_bridge_len:
            # check whether comp connects two (or more) distinct head atoms
            neighbor_heads=set()
            for i in comp:
                ai=mol.GetAtomWithIdx(i)
                for nb in ai.GetNeighbors():
                    if nb.GetIdx() in head_set:
                        neighbor_heads.add(nb.GetIdx())
            if len(neighbor_heads) >= 2:
                bridges.append(comp)
    if bridges:
        for comp in bridges:
            head_set |= comp
    return head_set, tails_list

# ------------- Benzene ring integrity (avoid fragmentation) ------------- #
def enforce_full_benzenes(mol, head_set, tails_list):
    """Ensure any 6-carbon aromatic ring is not split between head and tail.
    If any atom of the ring is in tail, move the entire ring to that tail (first intersecting tail)."""
    if not tails_list:
        return head_set, tails_list
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring)==6 and all(mol.GetAtomWithIdx(i).GetAtomicNum()==6 and mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            rset=set(ring)
            target=None
            for t in tails_list:
                if t & rset:
                    target=t; break
            if target is not None:
                if not rset.issubset(target):
                    target |= rset
                head_set -= rset
                # remove from other tails for consistency
                for t in tails_list:
                    if t is not target:
                        t -= rset
    return head_set, tails_list

# ------------- Promotion of carbohydrate-like rings ------------- #
def promote_carbohydrate_rings(mol, head_set, tails_list):
        # """Detect 5–6 member rings rich in oxygen (>=3 O in ring or as direct substituents) and move them entirely to head.
        # Includes:
        #     - Ring atoms
        #     - Direct neighboring oxygens (OH, anomeric)
        #     - CH2OH carbon (exocyclic carbon bound to a ring carbon and an oxygen) and its adjacent oxygen(s)
        # Does not include the alkyl chain carbon bound via anomeric oxygen (remains in tail), but includes the glycosidic oxygen.
        # """
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    head_set = set(head_set)
    # Precompute tail union to distinguish CH2OH vs long chain
    tail_union = set().union(*tails_list) if tails_list else set()
    for ring in rings:
        if len(ring) not in (5,6):
            continue
        ring_set = set(ring)
        # count oxygens in the ring
        oxy_in_ring = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum()==8)
        # adjacent (exocyclic) oxygens directly bound to ring atoms
        oxy_exo = set()
        for i in ring:
            a = mol.GetAtomWithIdx(i)
            for nb in a.GetNeighbors():
                if nb.GetIdx() not in ring_set and nb.GetAtomicNum()==8:
                    oxy_exo.add(nb.GetIdx())
        oxy_total = oxy_in_ring + len(oxy_exo)
        if oxy_total < 3:
            continue  # not clearly carbohydrate-like
        # Build sugar atom set
        sugar_atoms = set(ring_set) | oxy_exo
        # Add CH2OH carbons (carbon outside the ring connected to a ring carbon and an oxygen)
        for i in ring:
            a = mol.GetAtomWithIdx(i)
            for nb in a.GetNeighbors():
                nidx = nb.GetIdx()
                if nidx in ring_set:
                    continue
                if nb.GetAtomicNum()==6:
                    # has oxygen neighbor outside ring -> likely CH2OH
                    has_exo_O = any(o.GetAtomicNum()==8 and o.GetIdx() not in ring_set for o in nb.GetNeighbors())
                    if has_exo_O:
                        # ensure this is not the beginning of the long tail (heuristic: >2 consecutive linear carbons away from ring)
                        # count linear chain from nb excluding paths that return to ring
                        linear_len=0
                        visited={i}
                        stack=[(nb,0)]
                        while stack:
                            cur,depth=stack.pop()
                            if depth>8: break
                            linear_len=max(linear_len, depth)
                            cur_atom=mol.GetAtomWithIdx(cur.GetIdx())
                            for nn in cur_atom.GetNeighbors():
                                nid=nn.GetIdx()
                                if nid in ring_set or nid==i or nid in visited:
                                    continue
                                if nn.GetAtomicNum()==6 and hetero_neighbors(nn)==0:
                                    visited.add(nid)
                                    stack.append((nn,depth+1))
                        # If explored linear chain is very short (<3), assume CH2OH and not main tail
                        if linear_len < 3:
                            sugar_atoms.add(nidx)
                            for o in nb.GetNeighbors():
                                if o.GetAtomicNum()==8:
                                    sugar_atoms.add(o.GetIdx())
        # Update partition: remove sugar atoms from tails and add to head
        intersect_tail = any(sugar_atoms & t for t in tails_list)
        if intersect_tail:
            for t in tails_list:
                if sugar_atoms & t:
                    t -= sugar_atoms
        head_set |= sugar_atoms
    return head_set, tails_list

# ------------- Promotion of open polyol chains ------------- #
def promote_polyol_chains(mol, head_set, tails_list, min_chain_carbons=3):
        # """Promote/rescue open poly-hydroxylated chains (polyols) into head.
        # Relaxed version to avoid losing intermediate carbon without direct OH (e.g., HO-CH2-CH(-OH)-CH2-CH(-OH)-CH2-OH).

        # Strategy:
        #     1. Traverse aliphatic carbons (non-aromatic) not yet in head (may be in tail or free) and build components.
        #     2. For each C, mark as 'directly oxygenated' (>=1 O neighbor) or 'indirect' (no O but has C neighbor that has O).
        #     3. Acceptance criteria for a polyol component:
        #          - len(comp) >= min_chain_carbons
        #          - direct_oxygenated >= 2
        #          - (direct_oxygenated / len(comp)) >= 0.4  (minimum density)
        #          - (direct_oxygenated + indirect_oxygenated) == len(comp)  (all C within one step of an O)
        #          - Unique associated oxygens >= max(2, ceil(len(comp)*0.6))
        # If met, move carbons + associated oxygens to head, removing them from any tail.
        # """
    import math
    head_set = set(head_set)
    tail_sets = tails_list
    visited=set()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum()!=6 or atom.GetIsAromatic() or idx in head_set or idx in visited:
            continue
        # Build connected aliphatic-carbon component
        stack=[idx]; comp=set()
        while stack:
            i=stack.pop()
            if i in comp: continue
            ai=mol.GetAtomWithIdx(i)
            if ai.GetAtomicNum()!=6 or ai.GetIsAromatic():
                continue
            comp.add(i)
            for nb in ai.GetNeighbors():
                j=nb.GetIdx()
                if nb.GetAtomicNum()==6 and not nb.GetIsAromatic() and j not in comp and j not in head_set:
                    stack.append(j)
        visited |= comp
        if len(comp) < min_chain_carbons:
            continue
        # Oxygen analysis
        direct_oxy_c=set()
        oxy_atoms=set()
        for cidx in comp:
            ca=mol.GetAtomWithIdx(cidx)
            for nb in ca.GetNeighbors():
                if nb.GetAtomicNum()==8:
                    direct_oxy_c.add(cidx)
                    oxy_atoms.add(nb.GetIdx())
        if len(direct_oxy_c) < 2:
            continue
        indirect_c=set()
        for cidx in comp:
            if cidx in direct_oxy_c:
                continue
            ca=mol.GetAtomWithIdx(cidx)
            if any(nb.GetIdx() in direct_oxy_c for nb in ca.GetNeighbors() if nb.GetAtomicNum()==6):
                indirect_c.add(cidx)
        # Density and coverage
        if len(direct_oxy_c)/len(comp) < 0.4:
            continue
        if len(direct_oxy_c) + len(indirect_c) != len(comp):
            continue
        # Unique oxygen count density
        if len(oxy_atoms) < max(2, math.ceil(len(comp)*0.6)):
            continue
        # Promote
        to_move = comp | oxy_atoms
        # Remove from tails
        if any(to_move & t for t in tail_sets):
            for t in tail_sets:
                if to_move & t:
                    t -= to_move
        head_set |= to_move
    return head_set, tails_list

# ------------- Split head into functional fragments and cleanup ------------- #
def refine_heads_split(mol, head_set, tails_list):
        # """Generic version: split head into connected components avoiding special patterns.
        # Rules:
        #     - Compute connected components within head_set.
        #     - Aromatic bridging oxygen (O with two aromatic C neighbors) is moved to tail (not considered head).
        #     - Single-atom components that are isolated neutral oxygen are moved to tail.
        #     - Remaining components are returned as Head1, Head2, ...
        # """
    head_set = set(head_set)
    # 1. Move aromatic bridge oxygens to tail
    bridge_os=[]
    for idx in list(head_set):
        a=mol.GetAtomWithIdx(idx)
        if a.GetAtomicNum()==8:
            nbs=a.GetNeighbors()
            if len(nbs)==2 and all(nb.GetAtomicNum()==6 and nb.GetIsAromatic() for nb in nbs):
                bridge_os.append(idx)
    for oidx in bridge_os:
        a=mol.GetAtomWithIdx(oidx)
        neigh_ids={nb.GetIdx() for nb in a.GetNeighbors()}
        placed=False
        for t in tails_list:
            if t & neigh_ids:
                t.add(oidx); placed=True; break
        if not placed:
            if tails_list:
                tails_list[0].add(oidx)
            else:
                tails_list.append({oidx})
        head_set.discard(oidx)

    # 2. Connected components in head
    def components(atom_indices):
        comps=[]; visited=set(); atom_indices=set(atom_indices)
        for i in list(atom_indices):
            if i in visited: continue
            stack=[i]; comp=set()
            while stack:
                j=stack.pop()
                if j in comp: continue
                if j not in atom_indices: continue
                comp.add(j)
                aj=mol.GetAtomWithIdx(j)
                for nb in aj.GetNeighbors():
                    if nb.GetIdx() in atom_indices:
                        stack.append(nb.GetIdx())
            visited |= comp
            comps.append(comp)
        return comps
    comps = components(head_set)

    head_fragments=[]
    for comp in comps:
        if len(comp)==1:
            idx=next(iter(comp))
            a=mol.GetAtomWithIdx(idx)
            # Move isolated neutral O to tail
            # Avoid degrading if it would be the only existing head
            if a.GetAtomicNum()==8 and a.GetFormalCharge()==0 and not (len(head_set)==1):
                # try adding to a tail with a neighbor
                neigh_ids={nb.GetIdx() for nb in a.GetNeighbors()}
                placed=False
                for t in tails_list:
                    if t & neigh_ids:
                        t.add(idx); placed=True; break
                if not placed:
                    if tails_list:
                        tails_list[0].add(idx)
                    else:
                        tails_list.append({idx})
                head_set.discard(idx)
                continue  # do not keep as head fragment
        head_fragments.append(comp)

    # 3. Recompute head_set from valid fragments
    head_set = set().union(*head_fragments) if head_fragments else set()
    return head_set, head_fragments, tails_list

# ------------- Raw material inference ------------- #
def count_ethoxylate_units(smiles:str) -> int:
    # heuristic: count non-overlapping occurrences of OCC (approx -O-CH2-CH2-)
    return smiles.count('OCC')

def classify_tail_basic(mol, tail_atoms:set) -> dict:
    """
    New logic: TailChainLength = total number of carbons across all tails (sum),
    not the longest connected chain. For two tails of 4C and 4C => TailChainLength = 8.
    """
    if not tail_atoms:
        return {'TailChainLength':0,'TailUnsaturation':0,'TailAromatic':False}

    # Total carbons across all tails
    total_carbons = sum(1 for i in tail_atoms if mol.GetAtomWithIdx(i).GetAtomicNum()==6)

    # Unsaturations (count unique double/triple/aromatic bonds within tail)
    unsat_bonds = set()
    aromatic = False
    for idx in tail_atoms:
        a = mol.GetAtomWithIdx(idx)
        if a.GetIsAromatic():
            aromatic = True
        for b in a.GetBonds():
            other = b.GetOtherAtom(a)
            if other.GetIdx() in tail_atoms and b.GetBondType().name in ('DOUBLE','TRIPLE','AROMATIC'):
                pair = tuple(sorted((idx, other.GetIdx())))
                unsat_bonds.add(pair)
    unsat = len(unsat_bonds)

    # Carbons adjacent to any hetero (for fatty vs functionalized differentiation)
    hetero_adj = sum(
        1 for idx in tail_atoms
        if mol.GetAtomWithIdx(idx).GetAtomicNum()==6 and any(
            nb.GetAtomicNum() not in (1,6) for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
        )
    )

    return {
        'TailChainLength': total_carbons,
        'TailUnsaturation': unsat,
        'TailAromatic': aromatic
    }

def classify_head_basic(mol, head_atoms:set, head_smiles:str) -> dict:
    hetero=[i for i in head_atoms if mol.GetAtomWithIdx(i).GetAtomicNum() not in (1,6)]
    nO=sum(1 for i in hetero if mol.GetAtomWithIdx(i).GetAtomicNum()==8)
    nN=sum(1 for i in hetero if mol.GetAtomWithIdx(i).GetAtomicNum()==7)
    nS=sum(1 for i in hetero if mol.GetAtomWithIdx(i).GetAtomicNum()==16)
    etho_units = count_ethoxylate_units(head_smiles)
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        if set(ring).issubset(head_atoms):
            ring_atoms.update(ring)
    # sugar/polyol heuristic
    sugar_like = (nO >= MIN_SUGAR_OH and len(ring_atoms)>=5)
    polyol_like = (nO >= MIN_POLYOL_OH and etho_units==0 and not sugar_like)
    # Head charge-type classification (nonionic vs zwitterionic vs cationic/anionic)
    head_formal_charges = [mol.GetAtomWithIdx(i).GetFormalCharge() for i in head_atoms]
    any_pos = any(c>0 for c in head_formal_charges)
    any_neg = any(c<0 for c in head_formal_charges)
    net_charge = sum(head_formal_charges)
    if any_pos and any_neg and net_charge == 0:
        head_charge_type = 'zwitterionic'
    elif net_charge == 0 and not any_pos and not any_neg:
        head_charge_type = 'nonionic'
    elif net_charge > 0:
        head_charge_type = 'cationic'
    elif net_charge < 0:
        head_charge_type = 'anionic'
    else:
        # Residual case: mixed charges with net 0 that did not satisfy the first if (e.g., poorly assigned formal charges)
        head_charge_type = 'mixed'
    return {'HeadO':nO,'HeadN':nN,'HeadS':nS,'EthoxylateUnits':etho_units
            ,'Head_charge_type': head_charge_type}

DRAW_DIR = None  # set in main if requested

from rdkit.Chem import rdMolDescriptors

def _longest_consecutive_apolar_chain(mol, tail_atoms:set, mode:str='fast') -> int:
    """Length of the longest consecutive carbon chain within the tail.
    fast: uses subgraph diameter (2 BFS) -> O(N+E)
    full: fallback to exhaustive DFS for maximum precision (expensive).
    """
    carbon_tail = [i for i in tail_atoms if mol.GetAtomWithIdx(i).GetAtomicNum()==6]
    if not carbon_tail:
        return 0
    # Build neighbor list (tail carbons only)
    neigh = {i:[n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in tail_atoms and mol.GetAtomWithIdx(n.GetIdx()).GetAtomicNum()==6] for i in carbon_tail}
    if mode != 'full':  # FAST approximate diameter (exact for unweighted graphs)
        from collections import deque
        def bfs(start):
            dist={start:0}; q=deque([start])
            while q:
                v=q.popleft()
                for w in neigh[v]:
                    if w not in dist:
                        dist[w]=dist[v]+1; q.append(w)
            # returns farthest node and distances
            far=max(dist, key=lambda k: dist[k])
            return far, dist
        any_node = carbon_tail[0]
        far1,_ = bfs(any_node)
        far2, dist2 = bfs(far1)
        diameter_edges = max(dist2.values()) if dist2 else 0
        # chain length in number of carbons = edges + 1 (if at least one node)
        return diameter_edges + 1 if carbon_tail else 0
    # FULL (original exhaustive DFS)
    best=1
    termini=[i for i in carbon_tail if len(neigh[i])<=1]
    starts = termini if termini else carbon_tail
    stack=[(s, None, 1) for s in starts]
    seen_states=set()
    while stack:
        cur, prev, length = stack.pop()
        if length>best: best=length
        key=(cur, prev, length)
        if key in seen_states: continue
        seen_states.add(key)
        for nb in neigh[cur]:
            if nb==prev: continue
            stack.append((nb, cur, length+1))
    return best

def _tail_branching_index(mol, tail_atoms:set) -> float:
    carbon_tail = [i for i in tail_atoms if mol.GetAtomWithIdx(i).GetAtomicNum()==6]
    if not carbon_tail:
        return 0.0
    branch_points=0
    for i in carbon_tail:
        a=mol.GetAtomWithIdx(i)
        deg=sum(1 for n in a.GetNeighbors() if n.GetIdx() in carbon_tail)
        if deg>=3:
            branch_points+=1
    return branch_points/len(carbon_tail)

def _atomwise_tpsa(mol):
    # Returns list of per-atom TPSA contributions aligned with atom indices
    try:
        contribs, _ = rdMolDescriptors._CalcTPSAContribs(mol)  # type: ignore[attr-defined]
        # contribs is list of floats
        return contribs
    except Exception:
        return [0.0]*mol.GetNumAtoms()

def _crippen_atom_contribs(mol):
    try:
        contribs = rdMolDescriptors._CalcCrippenContribs(mol)  # type: ignore[attr-defined]
        # list of (logP, MR)
        return contribs
    except Exception:
        return [(0.0,0.0)]*mol.GetNumAtoms()

def compute_partition_descriptors(mol, head_atoms:set, tail_atoms:set, mode:str='fast') -> dict:
        """Reduced version: only base structural metrics requested for tail and head charge.
            Returns:
                - Tail_length (maximum consecutive carbon chain length)
                - Tail_branching_index (branch points / # carbon atoms in tail)
                - Head_charge (sum of formal charges in head)
            """
        head_atoms=set(head_atoms); tail_atoms=set(tail_atoms)
        tail_length=_longest_consecutive_apolar_chain(mol, tail_atoms, mode=mode)
        tail_branch=_tail_branching_index(mol, tail_atoms)
        head_charge=sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in head_atoms)
        return {
                'Tail_length': tail_length,
                'Tail_branching_index_calc': tail_branch,
                'Head_charge': head_charge
        }

# ---------------- Region feature set (TPSA/logP/charge range) using Head and TailCombined SMILES ---------------- #
def compute_region_features_from_smiles(head_smiles:str, tail_smiles:str) -> dict:
    """Compute TPSA, logP, and charge range for head and tail using complete fragment SMILES.
    head - tail for ΔTPSA; tail - head for ΔlogP; head - tail for ΔChargeRange.
    """
    from rdkit.Chem import Crippen
    hm = Chem.MolFromSmiles(head_smiles) if head_smiles else None
    tm = Chem.MolFromSmiles(tail_smiles) if tail_smiles else None
    # TPSA
    def safe_tpsa(m):
        if not m: return 0.0
        try: return rdMolDescriptors.CalcTPSA(m)  # type: ignore[attr-defined]
        except Exception: return 0.0
    # logP
    def safe_logp(m):
        if not m: return 0.0
        try: return Crippen.MolLogP(m)
        except Exception: return 0.0
    # Charge range (Gasteiger)
    def charge_range(m):
        if not m: return 0.0
        try:
            AllChem.ComputeGasteigerCharges(m)
            charges=[]
            for i in range(m.GetNumAtoms()):
                a=m.GetAtomWithIdx(i)
                try: q=a.GetDoubleProp('_GasteigerCharge')
                except Exception: q=0.0
                charges.append(q if q is not None else 0.0)
            return (max(charges)-min(charges)) if len(charges)>1 else 0.0
        except Exception:
            return 0.0
    head_tpsa = safe_tpsa(hm)
    tail_tpsa = safe_tpsa(tm)
    head_logp = safe_logp(hm)
    tail_logp = safe_logp(tm)
    head_crange = charge_range(hm)
    tail_crange = charge_range(tm)
    return {
        'TPSA_head': head_tpsa,
        'TPSA_tail': tail_tpsa,
        'Delta_TPSA': head_tpsa + tail_tpsa,
        'logP_head': head_logp,
        'logP_tail': tail_logp,
        'Delta_logP': head_logp +tail_logp ,
        'ChargeRange_head': head_crange,
        'ChargeRange_tail': tail_crange,
        'Delta_ChargeRange': head_crange + tail_crange
    }

# ------------- CPP (Critical Packing Parameter) related features ------------- #
def compute_cpp_features(head_smiles:str, tail_smiles:str, tail_chain_length:int) -> dict:
        # """Compute parameters for CPP with consistent units.

        # Units and formulas (Tanford / Evans):
        #     - Labute ASA (rdMolDescriptors.CalcLabuteASA) returns area in Å^2.
        #         a0_nm2 = ASA_Å2 * 0.01 (1 Å^2 = 0.01 nm^2)
        #     - Extended chain length: l_c(Å) = 1.5 + 1.265 * n_c  -> l_c_nm = l_c(Å)*0.1

        # """
    from rdkit.Chem import Crippen, Lipinski
    head_mol = Chem.MolFromSmiles(head_smiles) if head_smiles else None
    tail_mol = Chem.MolFromSmiles(tail_smiles) if tail_smiles else None

    # Head area (Å^2)
    def labute(m):
        if not m: return 0.0
        try:
            return rdMolDescriptors.CalcLabuteASA(m)  # type: ignore[attr-defined]
        except Exception:
            return 0.0
    head_area_A2 = labute(head_mol)

    # Count carbons in tail fragment 
    tail_total_carbons = 0
    if tail_mol:
        for a in tail_mol.GetAtoms():
            if a.GetAtomicNum()==6:
                tail_total_carbons += 1

    # Molar refractivity proxy (legacy)
    def molmr(m):
        if not m: return 0.0
        try: return Crippen.MolMR(m)
        except Exception: return 0.0
    tail_volume_mr = molmr(tail_mol)

    # Extended length (Å) and nm
    if tail_chain_length and tail_chain_length>0:
        tail_length_ext_A = 1.5 + 1.265 * tail_chain_length
    else:
        tail_length_ext_A = 0.0
    tail_length_ext_nm = tail_length_ext_A * 0.1 if tail_length_ext_A>0 else 0.0


    # Effective area in nm^2
    a0_nm2 = head_area_A2 * 0.01 if head_area_A2>0 else 0.0

    # HBA/HBD
    def hba(m):
        if not m: return 0
        try: return Lipinski.NumHAcceptors(m)
        except Exception: return 0
    def hbd(m):
        if not m: return 0
        try: return Lipinski.NumHDonors(m)
        except Exception: return 0

    return {
        # Legacy values (to avoid breaking existing pipelines)
        'Head_area': head_area_A2,                  # Å^2
        'Tail_volume': tail_volume_mr,              # MR proxy (not used in new CPP)
        'Tail_length_extended': tail_length_ext_A,  # Å
        'HBA_head': hba(head_mol),
        'HBD_head': hbd(head_mol),
        # New fields with explicit units
        'Head_area_nm2': a0_nm2,
        'Tail_nC_longest': tail_chain_length,
        'Tail_total_carbons': tail_total_carbons}

def draw_partition(mol, head_atoms, tail_atoms_list, out_path=None):
    """Generate the combo image (head vs tail side by side).

    If `out_path` is provided, saves PNG and returns the path.
    If `out_path` is None, returns a PIL Image object (or None if unavailable).
    """
    from rdkit.Chem import Draw as _Draw
    base = None
    if out_path:
        base, _ = os.path.splitext(out_path)
    head_atoms = set(head_atoms)
    tails = [set(t) for t in tail_atoms_list]
    union_tail = set().union(*tails) if tails else set()
    head_color=(0.12,0.35,0.9)  # blue
    tail_palette=[(0.95,0.55,0.1),(0.15,0.65,0.25),(0.7,0.3,0.75)]  # orange, green, purple

    # Copy for coordinates
    m = Chem.Mol(mol)
    try:
        from rdkit.Chem import rdDepictor
        rdDepictor.Compute2DCoords(m)
    except Exception:
        pass

    # Create individual images in memory (not saved) -----------------
    def make_head_img():
        if not head_atoms:
            return None
        img = _Draw.MolToImage(m, size=(400,400), highlightAtoms=list(head_atoms), highlightAtomColors={i:head_color for i in head_atoms})
        try:
            from PIL import ImageDraw
            d=ImageDraw.Draw(img); d.text((4,4), f"Head ({len(head_atoms)})", fill=(0,0,0))
        except Exception:
            pass
        return img

    def make_tail_img():
        if not union_tail:
            return None
        tail_color_map={}
        for ti, tset in enumerate(tails):
            col=tail_palette[ti % len(tail_palette)]
            for idx in tset:
                tail_color_map.setdefault(idx,col)
        img = _Draw.MolToImage(m, size=(400,400), highlightAtoms=list(union_tail), highlightAtomColors=tail_color_map)
        try:
            from PIL import ImageDraw
            d=ImageDraw.Draw(img); d.text((4,4), f"Tail ({len(union_tail)})", fill=(0,0,0))
        except Exception:
            pass
        return img

    head_img = make_head_img()
    tail_img = make_tail_img()

    # Build combo ---------------------------------------------------------
    try:
        from PIL import Image
        if head_img and tail_img:
            combo = Image.new('RGB', (head_img.width + tail_img.width, max(head_img.height, tail_img.height)), (255,255,255))
            combo.paste(head_img, (0,0))
            combo.paste(tail_img, (head_img.width,0))
        elif head_img:  # head only
            combo = head_img
        elif tail_img:  # tail(s) only
            combo = tail_img
        else:  # nothing to highlight; draw simple molecule
            combo = _Draw.MolToImage(m, size=(400,400))
        if base:
            combo_path = base + '_combo.png'
            combo.save(combo_path)
            return combo_path
        return combo
    except Exception:
        # Fallback without Pillow: try a single combined image (if RDKit colors suffice)
        try:
            highlight_atoms = list(head_atoms | union_tail)
            color_map = {i:head_color for i in head_atoms}
            # Assign tail colors
            for ti, tset in enumerate(tails):
                col=tail_palette[ti % len(tail_palette)]
                for idx in tset:
                    if idx not in color_map:
                        color_map[idx]=col
            img = _Draw.MolToImage(m, size=(500,500), highlightAtoms=highlight_atoms, highlightAtomColors=color_map)
            if base:
                combo_path = base + '_combo.png'
                img.save(combo_path)
                return combo_path
            return img
        except Exception as e:
            # Last resort: MolToFile without highlighting
            try:
                if base:
                    combo_path = base + '_combo.png'
                    _Draw.MolToFile(m, combo_path, size=(400,400))
                    return combo_path
                return _Draw.MolToImage(m, size=(400,400))
            except Exception:
                print('[WARN] Could not save combo image:', e)
                return None

def process_smiles(smiles, render_images=True):
    rows=[]
    for smi in smiles:
        mol = mol_ok(smi)
        if not mol: continue
        parts = split_head_tail(mol)
        
        tail_atoms = set().union(*parts['TailAtomsList']) if parts['TailAtomsList'] else set()
        tail_info = classify_tail_basic(mol, tail_atoms)
        head_info = classify_head_basic(mol, parts['HeadAtoms'], parts['Head'])

        part_desc = compute_partition_descriptors(mol, parts['HeadAtoms'], tail_atoms, mode=( 'fast'))
        # Region feature set (independent of part_desc naming for users)
        
        region_feats = compute_region_features_from_smiles(parts['Head'], parts['TailCombined'])
        # Chain length for Tanford: use TailChainLength (already calculated) or 0
        chain_len_for_cpp = tail_info.get('TailChainLength') or 0
        cpp_feats = compute_cpp_features(parts['Head'], parts['TailCombined'], chain_len_for_cpp)

        # Molecular weight of head and whole molecule
        try:
            mw_full = Descriptors.MolWt(mol)
        except Exception:
            mw_full = None
        try:
            head_mol = Chem.MolFromSmiles(parts['Head']) if parts['Head'] else None
            mw_head = Descriptors.MolWt(head_mol) if head_mol else None
        except Exception:
            mw_head = None
        row = {
            'Surfactant': smi,
            'CounterionList': ';'.join(parts.get('CounterionList', [])) if parts.get('CounterionList') else '',
            'Head': parts['Head'],
            'TailCombined': parts['TailCombined'],
            'HI': parts['HI'],
            'TailChainLength': tail_info['TailChainLength'],
            'TailUnsaturation': tail_info['TailUnsaturation'],
            'TailAromatic': tail_info['TailAromatic'],
            'HeadO': head_info['HeadO'],
            'HeadN': head_info['HeadN'],
            'HeadS': head_info['HeadS'],
            'EthoxylateUnits': head_info['EthoxylateUnits'],
            'Tail_branching_index': part_desc['Tail_branching_index_calc'],
            'Head_charge': part_desc['Head_charge'],
            'Head_charge_type': head_info.get('Head_charge_type'),

        }
        # Alias / user requested naming additions
        row.update({
            # Region metrics explicitly requested
            'TPSA_head': region_feats['TPSA_head'],
            'TPSA_tail': region_feats['TPSA_tail'],
            'ΔTPSA': region_feats['Delta_TPSA'],
            'logP_head': region_feats['logP_head'],
            'logP_tail': region_feats['logP_tail'],
            'ΔlogP': region_feats['Delta_logP'],
            'ChargeRange_head': region_feats['ChargeRange_head'],
            'ChargeRange_tail': region_feats['ChargeRange_tail'],
            'ΔChargeRange': region_feats['Delta_ChargeRange'],
            'Head_area': cpp_feats['Head_area'],
            'Tail_volume': cpp_feats['Tail_volume'],
            'Tail_length_extended': cpp_feats['Tail_length_extended'],
            'HBA_head': cpp_feats['HBA_head'],
            'HBD_head': cpp_feats['HBD_head'],

        })

        # Render images either to disk (if DRAW_DIR set) or in-memory (if render_images=True)
        if DRAW_DIR:
            import pathlib
            pathlib.Path(DRAW_DIR).mkdir(parents=True, exist_ok=True)
            img_name = f"mol_{len(rows)}.png"
            base_path = pathlib.Path(DRAW_DIR)/img_name
            _ = draw_partition(mol, parts['HeadAtoms'], parts['TailAtomsList'], str(base_path))
        elif render_images:
            img = draw_partition(mol, parts['HeadAtoms'], parts['TailAtomsList'], None)
            try:
                from PIL import Image
                if img is not None and hasattr(img, 'show'):
                    img.show(title=f"{smi}")
            except Exception:
                pass
        rows.append(row)
    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser(description='Head/Tail split & raw material inference for surfactants (heuristic).')
    parser.add_argument('--input', default=INPUT_EXCEL, help='Input Excel/CSV with surfactant SMILES (optional if --smiles is provided)')
    parser.add_argument('--smiles', nargs='*', help='One or more SMILES strings entered directly')
    parser.add_argument('--sheet', default=SHEET, help='Sheet name if Excel')
    parser.add_argument('--smiles-col', default=SMILES_COLUMN, help='Column containing SMILES strings')
    parser.add_argument('--output', default=None, help='Output Excel/CSV file') #default='output_excelfile.xlsx'
    parser.add_argument('--draw-dir', default=None, help='Directory to save colored PNGs (head/tail).')
    parser.add_argument('--descriptor-mode', choices=['none','fast','full'], default='fast', help='Partition descriptor level (none|fast|full). fast uses BFS; full uses exhaustive search for tail length.')
    parser.add_argument('--show-images', action='store_true', help='Render images in-memory instead of saving to disk')
    args = parser.parse_args()
    global DRAW_DIR

    # Resolve SMILES source: direct or file
    smiles = []
    if args.smiles:
        smiles = [s for s in args.smiles if mol_ok(s)]
    elif args.input:
        if not os.path.isfile(args.input):
            raise FileNotFoundError(args.input)
        if args.input.lower().endswith(('.xls','.xlsx')):
            df = pd.read_excel(args.input, sheet_name=args.sheet)
        else:
            df = pd.read_csv(args.input)
        if args.smiles_col not in df.columns:
            raise ValueError(f'Column {args.smiles_col} not in file')
        smiles = [s for s in df[args.smiles_col].astype(str) if mol_ok(s)]
    else:
        # Fallback: use built-in DEFAULT_SMILES so the script runs without console
        smiles = [s for s in DEFAULT_SMILES if mol_ok(s)]
    DRAW_DIR = args.draw_dir

    # If show-images is requested, do not force a draw directory
    if args.show_images:
        DRAW_DIR = None
        out = process_smiles(smiles, render_images=True)
    else:
        out = process_smiles(smiles, render_images=True)
    # Write or print output
    if args.output:
        if args.output.lower().endswith(('.xls', '.xlsx')):
            out.to_excel(args.output, index=False)
        else:
            out.to_csv(args.output, index=False)
        print('Output ->', args.output)
    else:
        # No output file specified; show results here
        try:
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', 120)
        except Exception:
            pass
        print('No output file specified; printing results below:')
        print(out.to_string(index=False))

    # Quick summary
    print('Processed', len(out), 'surfactants.')
    summary_cols = ['Surfactant']
    try:
        print(out[summary_cols].head())
    except Exception:
        print(out.head())

if __name__ == '__main__':
    main()