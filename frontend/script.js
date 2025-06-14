document.addEventListener('DOMContentLoaded', () => {
    // State
    let uploadedMolecules = []; // NEW: Store the actual data here
    let selectedPath = [];
    let allDiscoveredSmiles = new Set();
    let allExploredEdges = [];
    let currentParent = { Name: '', SMILES: '' };

    // DOM Elements
    const fileInput = document.getElementById('file-input');
    const fileNameDisplay = document.getElementById('file-name');
    const parentNameInput = document.getElementById('parent-name-input');
    const parentSmilesInput = document.getElementById('parent-smiles-input');
    const startAnalysisBtn = document.getElementById('start-analysis-btn');
    
    const branchingSection = document.getElementById('branching-section');
    const branchPointSelect = document.getElementById('branch-point-select');

    const candidatesGrid = document.getElementById('candidates-grid');
    const pathwayImageContainer = document.getElementById('pathway-image-container');
    const pathwayImage = document.getElementById('pathway-image');
    const pathwayPlaceholder = pathwayImageContainer.querySelector('.placeholder');
    const downloadBtn = document.getElementById('download-btn');
    
    const loadingDiv = document.getElementById('loading');
    const errorDiv = document.getElementById('error');

    // --- UI Update Functions ---
    const showLoading = (isLoading) => {
        loadingDiv.classList.toggle('hidden', !isLoading);
    };

    const showError = (message) => {
        errorDiv.textContent = message;
        errorDiv.classList.toggle('hidden', !message);
        if (message) {
            setTimeout(() => errorDiv.classList.add('hidden'), 5000); // Auto-hide after 5 seconds
        }
    };

    const updateStartButtonState = () => {
        startAnalysisBtn.disabled = !(uploadedMolecules.length > 0 && parentSmilesInput.value.trim());
    };

    const updateBranchingUI = () => {
        if (selectedPath.length > 1) {
            branchingSection.classList.remove('hidden');
            branchPointSelect.innerHTML = '';
            
            selectedPath.forEach((mol, index) => {
                const option = document.createElement('option');
                option.value = index;
                option.textContent = `${index + 1}: ${mol.Name}`;
                branchPointSelect.appendChild(option);
            });
            // Set the selected index to the last element
            branchPointSelect.selectedIndex = selectedPath.length - 1;
        } else {
            branchingSection.classList.add('hidden');
        }
    };

    // --- API Call Functions ---
    const handleApiCall = async (apiCall) => {
        showLoading(true);
        showError('');
        try {
            return await apiCall();
        } catch (err) {
            console.error(err);
            const errorMessage = err.response ? await err.response.json().then(d => d.detail) : 'An unexpected error occurred.';
            showError(`Error: ${errorMessage}`);
        } finally {
            showLoading(false);
        }
    };
    
    // --- Event Listeners ---
    fileInput.addEventListener('change', async () => {
        const file = fileInput.files[0];
        if (!file) return;

        const formData = new FormData();
        formData.append('file', file);
        
        await handleApiCall(async () => {
            const response = await fetch('/api/upload', { method: 'POST', body: formData });
            if (!response.ok) throw { response };
            
            uploadedMolecules = await response.json(); // Store the returned data
            
            fileNameDisplay.textContent = `File: ${file.name} (${uploadedMolecules.length} molecules loaded)`;
            updateStartButtonState();
        });
    });

    parentSmilesInput.addEventListener('input', updateStartButtonState);

    startAnalysisBtn.addEventListener('click', async () => {
        const initialParent = {
            Name: parentNameInput.value || 'Parent',
            SMILES: parentSmilesInput.value,
        };
        currentParent = initialParent;
        selectedPath = [initialParent];
        allDiscoveredSmiles.add(initialParent.SMILES);
        allExploredEdges = [];
        
        await handleApiCall(async () => {
            await fetchCandidatesAndPathway(initialParent);
        });
    });

    const fetchCandidatesAndPathway = async (parentMolecule) => {
        const excludedSmiles = Array.from(allDiscoveredSmiles);

        const analyzeRes = await fetch('/api/analyze', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ 
                parentSmiles: parentMolecule.SMILES, 
                selectedSmiles: excludedSmiles,
                molecules: uploadedMolecules // Send the full dataset
            })
        });
        if (!analyzeRes.ok) throw { response: analyzeRes };
        const analyzeData = await analyzeRes.json();
        renderCandidates(analyzeData.candidates);

        // Fetch pathway image using the new edge structure
        if (allExploredEdges.length > 0) {
            const pathwayRes = await fetch('/api/pathway-image', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ pathway_edges: allExploredEdges })
            });
            if (!pathwayRes.ok) throw { response: pathwayRes };
            const pathwayData = await pathwayRes.json();
            pathwayImage.src = pathwayData.image;
            pathwayImage.classList.remove('hidden');
            downloadBtn.classList.remove('hidden');
            if(pathwayPlaceholder) pathwayPlaceholder.classList.add('hidden');
        } else if (selectedPath.length === 1) {
            // If it's just the start, show only the parent image in the pathway view
             const imageRes = await fetch(`/api/molecule-image?smiles=${selectedPath[0].SMILES}`);
             if (!imageRes.ok) throw { response: imageRes };
             const imageData = await imageRes.json();
             pathwayImage.src = imageData.image;
             pathwayImage.classList.remove('hidden');
             downloadBtn.classList.remove('hidden');
             if(pathwayPlaceholder) pathwayPlaceholder.classList.add('hidden');
        }

        updateBranchingUI();
    };
    
    const renderCandidates = (candidates) => {
        candidatesGrid.innerHTML = ''; // Clear previous candidates
        if (candidates.length === 0) {
            const placeholder = document.createElement('p');
            placeholder.className = 'placeholder';
            placeholder.textContent = 'No further products found that meet the criteria.';
            candidatesGrid.appendChild(placeholder);
            return;
        }
        candidates.forEach(c => {
            const card = document.createElement('div');
            card.className = 'candidate-card';
            card.innerHTML = `
                <img src="${c.image}" alt="${c.Name}">
                <p>${c.Name}</p>
                <p>Similarity: ${c.ss.toFixed(3)}</p>
            `;
            card.addEventListener('click', () => selectCandidate(c));
            candidatesGrid.appendChild(card);
        });
    };

    const selectCandidate = async (candidate) => {
        const newEdge = { source: currentParent, target: candidate };
        allExploredEdges.push(newEdge);
        selectedPath.push(candidate);
        allDiscoveredSmiles.add(candidate.SMILES);
        const oldParent = currentParent;
        currentParent = candidate;

        await handleApiCall(async () => {
            // Update UI - only the top input boxes now, no 'Current Parent' section
            parentSmilesInput.value = candidate.SMILES;
            parentNameInput.value = candidate.Name;
            
            await fetchCandidatesAndPathway(currentParent);
        }).catch(() => {
            // If API fails, revert the state change
            allExploredEdges.pop();
            selectedPath.pop();
            allDiscoveredSmiles.delete(candidate.SMILES);
            currentParent = oldParent;
        });
    };

    const handleBranchSelection = async () => {
        const branchIndex = parseInt(branchPointSelect.value, 10);
        if (branchIndex === selectedPath.length - 1) {
            return;
        }

        const newPath = selectedPath.slice(0, branchIndex + 1);
        
        selectedPath = newPath;
        const branchParent = selectedPath[selectedPath.length - 1];
        currentParent = branchParent;
        
        await handleApiCall(async () => {
            // Update the UI inputs to reflect the new parent
            parentSmilesInput.value = branchParent.SMILES;
            parentNameInput.value = branchParent.Name;

            await fetchCandidatesAndPathway(branchParent);
        });
    }

    branchPointSelect.addEventListener('change', handleBranchSelection);

    downloadBtn.addEventListener('click', () => {
        const link = document.createElement('a');
        link.href = pathwayImage.src;
        link.download = 'reaction_pathway.png';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    });
}); 