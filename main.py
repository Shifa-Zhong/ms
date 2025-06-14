from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel, Field
import io
import base64
import os
from typing import List, Dict
from fastapi.staticfiles import StaticFiles
import tempfile

app = FastAPI()

# Mount the frontend directory
app.mount("/", StaticFiles(directory=os.path.join(os.path.dirname(os.path.realpath(__file__)), "frontend"), html=True), name="static")

# Helper functions from ms.py (to be moved here)
def count_atoms(smiles):
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if not molecule:
        return 0, 0, 0
    carbon_count = 0
    nitrogen_count = 0
    fluorine_count = 0
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_count += 1
        elif atom.GetSymbol() == 'N':
            nitrogen_count += 1
        elif atom.GetSymbol() == 'F':
            fluorine_count += 1
    return carbon_count, nitrogen_count, fluorine_count

def get_fp(smile):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np
    mol_mother = Chem.MolFromSmiles(smile)
    if not mol_mother:
        return np.array([])
    fp_mother = AllChem.GetHashedMorganFingerprint(mol_mother, 1, 2048)
    fp_mother = np.array(list(fp_mother)).astype('int32')
    return fp_mother

def tanimoto(arr1, arr2):
    import numpy as np
    if arr1.size == 0 or arr2.size == 0:
        return 0.0
    tanimoto_coefficient = np.dot(arr1, arr2) / (np.sum(arr1**2) + np.sum(arr2**2) - np.dot(arr1, arr2))
    return tanimoto_coefficient

def generate_molecule_image_base64(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    img = Draw.MolToImage(mol, size=(400, 400))
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    return f"data:image/png;base64,{img_str}"

class AnalyzeRequest(BaseModel):
    parentSmiles: str
    selectedSmiles: List[str]
    molecules: List[Dict]

class UploadResponse(BaseModel):
    sessionId: str
    fileName: str
    rowCount: int

class Molecule(BaseModel):
    Name: str
    SMILES: str
    ss: float
    image: str

class AnalyzeResponse(BaseModel):
    candidates: List[Molecule]

class PathwayMolecule(BaseModel):
    Name: str
    SMILES: str

class PathwayEdge(BaseModel):
    source: PathwayMolecule
    target: PathwayMolecule

class PathwayRequest(BaseModel):
    pathway_edges: List[PathwayEdge]

class ImageResponse(BaseModel):
    image: str

@app.post("/api/upload")
async def upload_file(file: UploadFile = File(...)):
    import pandas as pd
    try:
        file_type = file.filename.split('.')[-1]
        if file_type == 'xlsx':
            df = pd.read_excel(file.file)
        elif file_type == 'csv':
            df = pd.read_csv(file.file)
        else:
            raise HTTPException(status_code=400, detail="File format not supported.")
        
        if 'SMILES' not in df.columns or 'Name' not in df.columns:
            raise HTTPException(status_code=400, detail="File must contain 'SMILES' and 'Name' columns.")

        return df.to_dict(orient='records')
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/analyze", response_model=AnalyzeResponse)
async def analyze_molecules(request: AnalyzeRequest):
    import pandas as pd
    up_data = pd.DataFrame(request.molecules)
    
    carbon_m, nitrogen_m, fluorine_m = count_atoms(request.parentSmiles)
    fp_mother = get_fp(request.parentSmiles)

    candidates = []
    for i in up_data.index:
        smiles = up_data['SMILES'][i]
        if smiles not in request.selectedSmiles:
            carbon, nitrogen, fluorine = count_atoms(smiles)
            if carbon <= carbon_m and nitrogen <= nitrogen_m and fluorine <= fluorine_m:
                fp_son = get_fp(smiles)
                ss = tanimoto(fp_mother, fp_son)
                candidates.append({'SMILES': smiles, 'Name': up_data['Name'][i], 'ss': ss})
    
    sorted_candidates = sorted(candidates, key=lambda x: x['ss'], reverse=True)
    
    top_candidates = []
    for cand in sorted_candidates[:10]:
        image_b64 = generate_molecule_image_base64(cand['SMILES'])
        if image_b64:
            top_candidates.append(Molecule(Name=cand['Name'], SMILES=cand['SMILES'], ss=cand['ss'], image=image_b64))

    return AnalyzeResponse(candidates=top_candidates)

@app.get("/api/molecule-image", response_model=ImageResponse)
async def get_molecule_image(smiles: str):
    image_b64 = generate_molecule_image_base64(smiles)
    if not image_b64:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    return ImageResponse(image=image_b64)

@app.post("/api/pathway-image", response_model=ImageResponse)
async def generate_pathway_image(request: PathwayRequest):
    from rdkit import Chem
    from rdkit.Chem import Draw
    from PIL import Image
    import pydot
    if not request.pathway_edges:
        raise HTTPException(status_code=400, detail="Empty pathway provided.")

    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            pdot_graph = pydot.Dot(graph_type='digraph', rankdir='LR', splines='true', overlap='false', nodesep='0.8', ranksep='1.5')

            molecules_data = {}
            for edge in request.pathway_edges:
                for mol_data in [edge.source, edge.target]:
                    if mol_data.SMILES not in molecules_data:
                        mol = Chem.MolFromSmiles(mol_data.SMILES)
                        img = Draw.MolToImage(mol, size=(300, 300)) if mol else Image.new('RGBA', (300, 300))
                        
                        # Use a hash of the SMILES for a unique filename
                        img_filename = f"{hash(mol_data.SMILES)}.png"
                        img_path = os.path.join(temp_dir, img_filename)
                        img.save(img_path)

                        molecules_data[mol_data.SMILES] = {'name': mol_data.Name, 'image_path': img_path}

            # Add nodes to graph with images and labels
            for smiles, data in molecules_data.items():
                node = pydot.Node(
                    f'"{smiles}"',
                    label=data['name'],
                    shape='none', # Use image shape
                    image=data['image_path'],
                    labelloc='b', # Label at the bottom
                    fontname="Arial",
                    fontsize=24
                )
                pdot_graph.add_node(node)

            # Add edges
            for edge in request.pathway_edges:
                pdot_graph.add_edge(pydot.Edge(f'"{edge.source.SMILES}"', f'"{edge.target.SMILES}"'))

            # Create the final graph image using Graphviz
            png_bytes = pdot_graph.create_png()
            
            if not png_bytes:
                raise HTTPException(status_code=500, detail="Graphviz failed to generate PNG image.")

            img_str = base64.b64encode(png_bytes).decode()
            return ImageResponse(image=f"data:image/png;base64,{img_str}")

        except Exception as e:
            import traceback
            traceback.print_exc()
            raise HTTPException(status_code=500, detail=f"Failed to generate pathway graph image: {e}") 