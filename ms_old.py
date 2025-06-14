import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd
import numpy as np
from PIL import Image,ImageDraw, ImageFont
from rdkit.Chem import rdChemReactions
import networkx as nx
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="Welcome to MS-analysis",
    page_icon="💧",
    layout="wide",
    initial_sidebar_state="auto")
if 'rerun' not in st.session_state:
        st.session_state.rerun = 0
        
st.session_state.rerun+=1
st.write(st.session_state.rerun)
        
        
@st.cache_data(ttl=24*3600)
def load_data(file_name):
    file_type = str(file_name).split('.')[1]
    if file_type == 'xlsx':
        data = pd.read_excel(file_name)
    elif file_type == 'csv':
        data = pd.read_csv(file_name)
    else:
        st.warning('File format is not supported. Please upload the CSV or Excel file')
        st.stop()
    return data

def convert_df(df):
    return df.to_csv().encode('utf-8')

def count_atoms(smiles):
    molecule = Chem.MolFromSmiles(smiles)
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
    mol_mother = Chem.MolFromSmiles(smile)
    fp_mother = AllChem.GetHashedMorganFingerprint(mol_mother, 1, 2048)
    fp_mother = np.array(list(fp_mother)).astype('int32')
    return fp_mother

def tanimoto(arr1, arr2):
    tanimoto_coefficient = np.dot(arr1, arr2) / (np.sum(arr1**2) + np.sum(arr2**2) - np.dot(arr1, arr2))
    return tanimoto_coefficient

def update_new_data():
    sss = []
    ind = []
    carbon_m, nitrogen_m, fluorine_m = count_atoms(st.session_state.mother_smile)
    fp_mother = get_fp(st.session_state.mother_smile)
    for i in up_data.index:
        if up_data['SMILES'][i] not in st.session_state.selected_s:
            carbon, nitrogen, fluorine = count_atoms(up_data['SMILES'][i])
            if carbon <= carbon_m and nitrogen <= nitrogen_m and fluorine <= fluorine_m:
                ind.append(i)
                son = Chem.MolFromSmiles(up_data['SMILES'][i])
                fp_son = AllChem.GetHashedMorganFingerprint(son, 1, 2048)
                fp_son = np.array(list(fp_son)).astype('int32')
                ss = tanimoto(fp_mother, fp_son)
                sss.append(ss)
    new_data = up_data.loc[ind]
    new_data['ss'] = sss
    new_data.sort_values(by=['ss'], ascending=False, inplace=True)
    new_data.reset_index(drop=True, inplace=True)
    st.session_state.new_data = new_data
    
def update_product_image():
        if len(st.session_state.new_data)>=5:            
            st.session_state.options = [f"{st.session_state.new_data['Name'][i]}" for i in range(5)]
            for i in range(5):
                mol = Chem.MolFromSmiles(st.session_state.new_data['SMILES'][i])
                Draw.MolToImageFile(mol, f'Product_{i+1}.png', size=(400, 400))
        elif len(st.session_state.new_data)==4:
            st.session_state.options = [f"{st.session_state.new_data['Name'][i]}" for i in range(4)]
            for i in range(4):
                mol = Chem.MolFromSmiles(st.session_state.new_data['SMILES'][i])
                Draw.MolToImageFile(mol, f'Product_{i+1}.png', size=(400, 400))       
        elif len(st.session_state.new_data)==3:
            st.session_state.options = [f"{st.session_state.new_data['Name'][i]}" for i in range(3)]
            for i in range(3):
                mol = Chem.MolFromSmiles(st.session_state.new_data['SMILES'][i])
                Draw.MolToImageFile(mol, f'Product_{i+1}.png', size=(400, 400))               
        elif len(st.session_state.new_data)==2:
            st.session_state.options = [f"{st.session_state.new_data['Name'][i]}" for i in range(2)]
            for i in range(2):
                mol = Chem.MolFromSmiles(st.session_state.new_data['SMILES'][i])
                Draw.MolToImageFile(mol, f'Product_{i+1}.png', size=(400, 400))            
        elif len(st.session_state.new_data)==1:
            st.session_state.candidate_product_number=1
            st.session_state.options = [f"{st.session_state.new_data['Name'][0]}"]
            mol = Chem.MolFromSmiles(st.session_state.new_data['SMILES'][0])
            Draw.MolToImageFile(mol, f'Product_{1}.png', size=(400, 400))            
        else:
            st.session_state.options = [None]
            
def update_candiate_num():
    if len(st.session_state.new_data)>=5:
        st.session_state.candidate_product_number=5
    elif len(st.session_state.new_data)==4:
        st.session_state.candidate_product_number=4
    elif len(st.session_state.new_data)==3:
        st.session_state.candidate_product_number=3
    elif len(st.session_state.new_data)==2:
        st.session_state.candidate_product_number=2
    else:
        st.session_state.candidate_product_number=0
        
    
    
def updata_reaction_image(mother_name):    
    font = ImageFont.truetype("arial.ttf",15) 
    reaction_images = []
    r1 = st.session_state.selected_s[0]
    r2 = st.session_state.selected_s[1]
    r1_mol = Chem.MolFromSmiles(r1)
    r2_mol = Chem.MolFromSmiles(r2)
    smarts_step1 = Chem.MolToSmarts(r1_mol) + ">>" + Chem.MolToSmarts(r2_mol)
    reaction1 = rdChemReactions.ReactionFromSmarts(smarts_step1)
    img1 = Draw.ReactionToImage(reaction1)
    img1.convert("RGBA")
    width, height = img1.size
    new_height = height + 30  # Add extra space for names
    img_with_names = Image.new("RGBA", (width, new_height), (255, 255, 255, 255))

    # Paste the reaction image onto the new image
    img_with_names.paste(img1, (0, 0))

    # Add names under the reactant and product
    draw = ImageDraw.Draw(img_with_names)

    # Calculate positions for the names
    reactant_position = (width // 4, height + 5)  # Centered under the reactant
    product_position = (3 * width // 4, height + 5)  # Centered under the product
    name1 = up_data[up_data['SMILES']==r2]['Name'].values[0]

    # Draw the names
    draw.text(reactant_position, mother_name, font=font, fill=(0, 0, 0))
    draw.text(product_position, name1, font=font, fill=(0, 0, 0))
    reaction_images.append(img_with_names)

    if len(st.session_state.selected_s)>2:
        for sim in st.session_state.selected_s[2:]:
            name = up_data[up_data['SMILES']==sim]['Name'].values[0]
            p1_mol = Chem.MolFromSmiles(sim)
            smarts_step2 = ">>" + Chem.MolToSmarts(p1_mol)
            reaction2 = rdChemReactions.ReactionFromSmarts(smarts_step2)
            img2 = Draw.ReactionToImage(reaction2)
            img2 = img2.convert("RGBA")
            width, height = img2.size
            new_height = height + 30  # Add extra space for names
            img_with_names = Image.new("RGBA", (width, new_height), (255, 255, 255, 255))
            img_with_names.paste(img2, (0, 0))
            draw = ImageDraw.Draw(img_with_names)

            product_position = (2 * width // 4, height + 5)  # Centered under the product
            # Draw the names

            draw.text(product_position, f"{name}", font=font, fill=(0, 0, 0))

            reaction_images.append(img_with_names)

    total_width = sum(img.width for img in reaction_images)
    max_height = max(img.height for img in reaction_images)
    combined_img = Image.new('RGBA', (total_width, max_height))

    current_width = 0
    for img in reaction_images:
        combined_img.paste(img, (current_width, 0))
        current_width += img.width

    # Optionally save the combined image
    combined_img.save("sequential_reactions_from_list.png")


def updata_reaction_graph(mother_name):
    # Step 1: Create the graph
    G = nx.DiGraph()  # Use DiGraph for directed graph with arrows
    
    for i in range(len(st.session_state.selected_s)):
        mol = Chem.MolFromSmiles(st.session_state.selected_s[i])
        Draw.MolToImageFile(mol, f'PProduct_{i}.png', size=(400, 400))

    # Step 2: Add nodes with images and names
    image_paths = [f'PProduct_{i}.png' for i in range(len(st.session_state.selected_s))]
    node_names = [f'PProduct_{i}' for i in range(len(st.session_state.selected_s))]  # Names corresponding to each image

    for i, (image_path, name) in enumerate(zip(image_paths, node_names)):
        G.add_node(i, image=Image.open(image_path), name=name)

    # Step 3: Add directed edges between nodes
    # Step 3: Add directed edges between nodes for n nodes
    
    for i in range(len(st.session_state.selected_s) - 1):  # Loop through nodes and connect each node to the next
        G.add_edge(i, i + 1)

    # Step 4: Define the layout
    pos = nx.spring_layout(G)

    # Step 5: Draw the graph
    fig, ax = plt.subplots()

    # Define the desired image size on the plot (in axis coordinates)
    image_size = 0.5  # Adjust this value to increase or decrease image size

    # Draw the nodes with images and add labels
    for node in G.nodes:
        img = G.nodes[node]['image']

        # Resize the image
        img = img.resize((int(50/image_size), int(50/image_size)), Image.LANCZOS)

        x, y = pos[node]

        # Convert the PIL image to a numpy array and plot it with adjusted size
        img_array = np.array(img)
        ax.imshow(img_array, extent=(x-image_size/2, x+image_size/2, y-image_size/2, y+image_size/2), zorder=1)

        # Add text label to each node, positioned slightly below the image
        ax.text(x, y - image_size/2 - 0.1, G.nodes[node]['name'], horizontalalignment='center', fontsize=12, zorder=2)

    # Draw the edges with arrows
    nx.draw_networkx_edges(G, pos, ax=ax, arrowstyle='-|>', arrowsize=20, edge_color='black', width=2, connectionstyle='arc3,rad=0.1')

    # Hide the axis
    plt.axis('off')
    plt.savefig("sequential_reactions_from_list.png")
try:
    G
except NameError:
    G = nx.DiGraph()
    
def add_nodes_and_connect_sequential(node_list):
    global G
    G = st.session_state.G
    G.add_nodes_from(node_list)
    for i in range(len(node_list) - 1):
        G.add_edge(node_list[i], node_list[i + 1])        

   
    
def updata_reaction_graph():
    global G
    G = st.session_state.G
    node_colors = "grey"  # All nodes will be grey
    # Draw the graph
    pos = nx.spring_layout(G)  # You can change the layout for a different appearance
    nx.draw(G, pos, with_labels=True, node_color=node_colors, edge_color="black", node_size=1000, font_size=7, font_color="black", font_weight="bold", arrows=True)
    plt.axis('off')
    plt.savefig("sequential_reactions_from_list.png")


    
st.subheader('Please upload your files containing SMILES of possible products')
upload_file = st.file_uploader("Choose a csv or excel file containing molecules")
sample_data = load_data('data/example.xlsx')
if st.checkbox('Show an example of csv or excel file'):
    st.write(sample_data)
if upload_file is not None:
    file_type = str(upload_file.name).split('.')[1]
    if file_type == 'xlsx':
        up_data = pd.read_excel(upload_file)
    elif file_type == 'csv':
        up_data = pd.read_csv(upload_file)
    else:
        st.warning('File format is not supported. Please upload the CSV or Excel file')
        st.stop()
    with st.expander("Show your data"):
        st.write(up_data)
         
        
if upload_file is None:
    st.warning('You should first upload your files')
    st.stop()
else:   
    if 'selected_s' not in st.session_state:
        st.session_state.selected_s = []
    if 'mother_smile' not in st.session_state:
        st.session_state.mother_smile = []
    if 'mother_name' not in st.session_state:
        st.session_state.mother_name = []
    if 'path' not in st.session_state:
        st.session_state.path = 1
    if 'need_another_path' not in st.session_state:
        st.session_state.need_another_path = False
    if 'pathway_image' not in st.session_state:
        st.session_state.pathway_image = None
    if 'options' not in st.session_state:
        st.session_state.options=None
    if 'candidate_product_number' not in st.session_state:
        st.session_state.candidate_product_number=None
    if 'G' not in st.session_state:
        st.session_state.G = nx.DiGraph()

        

mother_s = st.text_input("What is the SMILES of parent chemical")
if mother_s:
    mother_name_user = st.text_input("How to name your parent chemical")

if not st.session_state.mother_smile:
    st.session_state.mother_smile = mother_s
    st.session_state.mother_name = mother_name_user
    st.session_state.selected_s.append(mother_s)  # Append the initial SMILES
else:
    # Use the existing values from session state
    mother_s = st.session_state.mother_smile
    mother_name = st.session_state.mother_name

try:
    carbon_m, nitrogen_m, fluorine_m = count_atoms(st.session_state.mother_smile)
except:
    st.warning('Please Check the SMILES of parent chemical')
    st.stop()

    
    
d1, d2 = st.columns(2)

with d1:
    col41, col42, col43, col44, col45, col46 = st.columns([1, 1, 1, 1,1,1])
    image_placeholder1 = col41.empty()
    image_placeholder2 = col42.empty()
    image_placeholder3 = col43.empty()
    image_placeholder4 = col44.empty()
    image_placeholder5 = col45.empty()
    image_placeholder6 = col46.empty()
    if st.session_state.mother_smile:
        fp_mother = get_fp(st.session_state.mother_smile)
        mol_mother = Chem.MolFromSmiles(st.session_state.mother_smile)
        Draw.MolToImageFile(mol_mother, 'mol_mother.png', size=(400, 400))
        image_placeholder1.image('mol_mother.png', caption=f'The image of {st.session_state.mother_name}')
        
        update_new_data()
        
        update_product_image()
        update_candiate_num()
        with st.expander("Show new data"):
            st.write(st.session_state.new_data)
        
        if st.session_state.candidate_product_number==5:
            image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
            image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
            image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
            image_placeholder5.image("Product_4.png", caption=f"{st.session_state.new_data['Name'][3]}")
            image_placeholder6.image("Product_5.png", caption=f"{st.session_state.new_data['Name'][4]}")
        
        if st.session_state.candidate_product_number==4:
            image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
            image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
            image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
            image_placeholder5.image("Product_4.png", caption=f"{st.session_state.new_data['Name'][3]}")
            image_placeholder6 = col46.empty()
        
        if st.session_state.candidate_product_number==3:
            image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
            image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
            image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
            image_placeholder5 = col45.empty()
            image_placeholder6 = col46.empty()
            
            
        if st.session_state.candidate_product_number==2:
            image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
            image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
            image_placeholder4 = col44.empty()
            image_placeholder5 = col45.empty()
            image_placeholder6 = col46.empty()
        if st.session_state.candidate_product_number==1:        
            image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
            image_placeholder3 = col43.empty()
            image_placeholder4 = col44.empty()
            image_placeholder5 = col45.empty()
            image_placeholder6 = col46.empty()
        if st.session_state.candidate_product_number==0:
            st.warning("Finished")
            
    
        # Dynamically update the options based on the new data
        
        selected_option = st.selectbox("What is your choice", options=st.session_state.options, index=None, placeholder="Select the degradation product")
        
        if selected_option is not None:
            st.session_state.mother_name = selected_option ##update mother name
            ind = st.session_state.new_data[st.session_state.new_data['Name'] == selected_option].index[0]
            st.session_state.mother_smile = st.session_state.new_data.loc[ind, 'SMILES'] ## update mother smile
            st.session_state.selected_s.append(st.session_state.mother_smile) ## update selected smile
            mol_mother = Chem.MolFromSmiles(st.session_state.mother_smile)
            Draw.MolToImageFile(mol_mother, 'mol_mother.png', size=(400, 400)) ## update the mother image
            image_placeholder1.image('mol_mother.png', caption=f'The image of {selected_option}')
    
            update_new_data() ## based on the new mother smile, updata the new data
    
            
    
            
            update_product_image()
            update_candiate_num()  
    
            if st.session_state.candidate_product_number==5:
                image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
                image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
                image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
                image_placeholder5.image("Product_4.png", caption=f"{st.session_state.new_data['Name'][3]}")
                image_placeholder6.image("Product_5.png", caption=f"{st.session_state.new_data['Name'][4]}")
        
            if st.session_state.candidate_product_number==4:
                image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
                image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
                image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
                image_placeholder5.image("Product_4.png", caption=f"{st.session_state.new_data['Name'][3]}")
                image_placeholder6 = col46.empty()
    
            if st.session_state.candidate_product_number==3:
                image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
                image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
                image_placeholder4.image("Product_3.png", caption=f"{st.session_state.new_data['Name'][2]}")
                image_placeholder5 = col45.empty()
                image_placeholder6 = col46.empty()
    
    
            if st.session_state.candidate_product_number==2:
                image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
                image_placeholder3.image("Product_2.png", caption=f"{st.session_state.new_data['Name'][1]}")
                image_placeholder4 = col44.empty()
                image_placeholder5 = col45.empty()
                image_placeholder6 = col46.empty()
            if st.session_state.candidate_product_number==1:        
                image_placeholder2.image("Product_1.png", caption=f"{st.session_state.new_data['Name'][0]}")
                image_placeholder3 = col43.empty()
                image_placeholder4 = col44.empty()
                image_placeholder5 = col45.empty()
                image_placeholder6 = col46.empty()
            if st.session_state.candidate_product_number==0:
                st.warning("Finished")
                
    
            selected_option = st.selectbox("What is your choice", st.session_state.options, index=None, placeholder="Select the degradation product")
    
            #updata_reaction_image(mother_name_user)
            
            nodes=[]
            for sim in st.session_state.selected_s:
                name = up_data[up_data['SMILES']==sim]['Name'].values[0]
                nodes.append(name)
            
            add_nodes_and_connect_sequential(nodes)
            
            updata_reaction_graph()
            
        
        
with d2:
    if st.session_state.mother_smile:
        if st.button("I need another pathway"):
            st.session_state.need_another_path = True
            
        if st.session_state.need_another_path:
            if len(st.session_state.selected_s) > 0:
                name_s=[]
                for sim in st.session_state.selected_s:
                    name = up_data[up_data['SMILES']==sim]['Name'].values[0]
                    name_s.append(name)
                new_start_molecule_name = st.selectbox("Which chemial do you want to start?", name_s, index=None, placeholder="Which position to start")
                if new_start_molecule_name:
                    ms = up_data[up_data['Name']==new_start_molecule_name]['SMILES'].values[0] ## located the selected smiles
                    ind = st.session_state.selected_s.index(ms)
                    st.session_state.selected_s = st.session_state.selected_s[0:ind+1] ### select all the smile from 1st to the selected one
                    
                    #updata_reaction_image(mother_name_user)## based on the updated selected_s to update the reaction image
                    update_candiate_num()
                    
                    #nodes=[]
                    #for sim in st.session_state.selected_s:
                     #   name = up_data[up_data['SMILES']==sim]['Name'].values[0]
                     #   nodes.append(name)
            
                    #add_nodes_and_connect_sequential(nodes)
                
                    #updata_reaction_graph()
                    
                    st.session_state.mother_smile = ms
                    st.session_state.mother_name = new_start_molecule_name
                    st.session_state.path += 1
                    st.rerun()                
                
            else:
                st.warning("You need at least one degradation product selected to start another pathway.")
                
        col51,_= st.columns([1,0.1])
        image_placeholder5 = col51.empty()
        try:
            image_placeholder5.image("sequential_reactions_from_list.png", caption=f"Pathway_{st.session_state.path}", output_format ="PNG", use_column_width="always")
        except:
            pass


    
    
            

                
        
            
            
                
                
                
                
