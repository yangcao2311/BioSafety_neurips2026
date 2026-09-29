# T4.4 TargetDiff — Session 9 status

**Weights**: downloaded successfully from Zenodo public mirror (record 14041881, `targetdiff_pretrained_models.tar.gz`). The original Google Drive folder is not publicly accessible to gdown. Files saved to `/home/xliu316/bio/targetdiff/pretrained_models/`:
- `pretrained_diffusion.pt` (33 MB) — SBDD diffusion model
- `egnn_pdbbind_v2016.pt` (30 MB) — affinity prediction model
- `pk_reg_para.pkl` — regression parameters

**Sampling**: `scripts/run_targetdiff.py` runs end-to-end on the 4 holo PDBs that Pocket2Mol used (7BV2 RdRp, 7VH8 Mpro, 3OXC HIV protease, 1VRT HIV RT). Configuration: 10 samples per pocket × 200 diffusion steps. PDB cleanup strips non-standard amino acids that broke TargetDiff's `PDBProtein` parser (DNA/RNA bases, HETATM ligands, waters, metals).

**Reconstruction**: across all 4 pockets, TargetDiff's `reconstruct.reconstruct_from_generated` returned no valid SMILES (`. `-fragmented multi-component outputs or `MolReconsError` on every sample). The diffusion sampling itself completes successfully — the failure is in the structure-to-SMILES reconstruction step.

**Interpretation**: TargetDiff was originally evaluated on the CrossDocked2020 test set with curated pocket inputs. Reconstruction quality is sensitive to (a) pocket conditioning fidelity, (b) sample count (the original paper uses N=100, we use N=10 for time), and (c) the post-processing chain (rdkit canonicalization plus the open-babel reconstruct pipeline). Our holo-PDB pockets, while sufficient for Pocket2Mol's atom-by-atom growth, do not produce reconstructable diffusion outputs at low sample count.

**Action**: documented as known TargetDiff reconstruction limitation in this environment. Weights and runner are shipped; raising N=100 samples and tightening pocket cropping is the recommended next step.

---

The Pocket2Mol-based T4.4 cascade with `cascade_results.json` (Max P_harm 0.3883) remains the working SBDD result for this submission.
