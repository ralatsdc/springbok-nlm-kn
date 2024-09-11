# NCBI Cell Pilot
## Motivation

The NCBI Cell pilot aims to create a comprehensive cell phenotype knowledgebase. This knowledgebase will capture details about cell phenotypes, making it available for research and integrating it with data about diseases and drugs. The goal is to discover new biomarkers and therapeutic targets.

We will use two main sources of knowledge:

1. **Single Cell Genomics Data**: Analyzing data from repositories like CZI CELLxGENE and NeMO to identify cell type-specific marker genes and differential expression patterns.

2. **Peer-Reviewed Publications**: Extracting information from PubMed Central (PMC) about single cell genomics experiments.

The extracted information will be linked with experiment metadata, such as species, sample sources, disease states, and responses to treatments. This data will be structured using semantic web technologies like RDF and OWL, and stored in graph databases like ArangoDB.

## Goals

This README provides guidance on using Python to:

0. **Set Up the Environment**: Prepare your development environment to properly run Jupyter Notebooks. This includes configuring environment variables and installing necessary Python dependencies. (See [Development Environment](#development-environment))

1. **Retrieve Human Lung Cell Data**: Use CELLxGENE to get dataset metadata and files. (See [Chapter-01-CELLxGENE.ipynb](ncbi-cell/ipynb/Chapter-01-CELLxGENE.ipynb))

2. **Identify Citations in PubMed**: Use NCBI E-Utilities to find PubMed citations for CELLxGENE datasets. (See [Chapter-02-E-Utilities.ipynb](ncbi-cell/ipynb/Chapter-02-E-Utilities.ipynb))

3. **Discover Marker Genes**: Use the NSForest package to find marker gene combinations from single-cell RNA sequencing data. (See [Chapter-03-NSforest.ipynb](ncbi-cell/ipynb/Chapter-03-NS-Forest.ipynb))

4. **Extract Structured Information**: Use OntoGPT to extract and structure information from text using large language models (LLMs) and ontology-based grounding. (See [Chapter-04-ontoGPT.ipynb](ncbi-cell/ipynb/Chapter-04-OntoGPT.ipynb))

5. **Visualize Data with ArangoDB**: Create, populate, and display results in a graph database. (See [Chapter-05-ArangoDB.ipynb](ncbi-cell/ipynb/Chapter-05-ArangoDB.ipynb))

## Development Environment

### Environment Variables

To use PubMed and OntoGPT, you need access keys called API keys. 

Here’s how to set up your API keys:

1. **Create a `.zshenv` File**

   You need to create a file named `.zshenv` in your home directory. This file will store your API keys. You can do this by opening a terminal and running the following command:

   ```sh
   touch ~/.zshenv
   ```

2. **Add Your API Keys**

   Open the `.zshenv` file using a text editor. Below is an example of how to use `nano`, but feel free to use any editor you prefer. 

   ```sh
   nano ~/.zshenv
   ```

   Add the following lines to the file. Replace `<YOUR_BIOPORTAL_API_KEY>`, `<YOUR_OPENAI_API_KEY>`, `<YOUR_NCBI_EMAIL>`, and `<YOUR_NCBI_API_KEY>` with your actual API keys and email:

   ```sh
   export BIOPORTAL_API_KEY=<YOUR_BIOPORTAL_API_KEY>
   export OPENAI_API_KEY=<YOUR_OPENAI_API_KEY>
   export NCBI_EMAIL=<YOUR_NCBI_EMAIL>
   export NCBI_API_KEY=<YOUR_NCBI_API_KEY>
   ```

   Save the file and exit the editor. In `nano`, you can do this by pressing `CTRL + X`, then `Y` to confirm, and `Enter` to save.


   See:

   - [API keys for E-Utilities](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
   - [OpenAI platform](https://openai.com/api/)
   - [OpenAI project API keys](https://platform.openai.com/api-keys)
   - [BioPortal](https://bioportal.bioontology.org/login)

3. **Apply the Changes**

   For the changes to take effect, you need to reload the `.zshenv` file. You can do this by running:

   ```sh
   source ~/.zshenv
   ```

   This command updates your terminal with the new environment variables you added.

**Important:** Do not share or commit your `.zshenv` file to this or any repositories. This file contains sensitive information.

### Python Dependencies

We use Poetry to manage dependencies. Before installing dependencies, you should have Python version 3.10, 3.11, or 3.12 installed. We choose to install 3.12 here. Follow these steps to set up your environment:

#### 1. Install Python 3.12

**For Windows:**

1. Download the Python 3.12 installer from the [official Python website](https://www.python.org/downloads/release/python-3125/).
2. Run the installer and ensure the option "Add Python 3.12 to PATH" is checked.
3. Click "Install Now" and follow the prompts to complete the installation.

**For macOS:**

1. Open Terminal.
2. Install Python 3.12 using Homebrew by running:
    ```sh
    brew install python@3.12
    ```
3. Verify the installation by running:
    ```sh
    python3.12 --version
    ```

If you see the error `brew: command not found`, it means Homebrew is not installed. To fix this:

1. **Install Homebrew**:
   Open Terminal and run:
   ```sh
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

2. **Verify Installation**:
   Check that Homebrew is installed by running:
   ```sh
   brew --version
   ```

   This should display the Homebrew version, confirming that it’s correctly installed. Run `brew install python@3.12` again to finish Python installation.

**For Linux:**

1. Open Terminal.
2. Install Python 3.12 using the package manager. For example, on Ubuntu, run:
    ```sh
    sudo apt update
    sudo apt install python3.12
    ```
3. Verify the installation by running:
    ```sh
    python3.12 --version
    ```

#### 2. Install Poetry

Poetry is used to manage project dependencies. Follow these steps to install Poetry:

1. Install Poetry in a virtual environment:

    ```sh
    python3.12 -m venv .poetry
    source .poetry/bin/activate
    python -m pip install -r .poetry.txt
    deactivate
    ```
   
Note that if you do not have Python 3.12 installed, you may need to modify these commands to point to the version of Python that is installed: e.g. `python3.10 -m venv .poetry`

#### 3. Install Project Dependencies

1. Create and activate a new virtual environment for the project:

    ```sh
    python3.12 -m venv .venv
    source .venv/bin/activate
    ```

2. Use Poetry to install the project's dependencies:

    ```sh
    .poetry/bin/poetry install
    ```

These steps will ensure you have Python 3.12 and all necessary dependencies installed to work with the Jupyter Notebooks in this project.


## Running the Notebooks

1. Activate your Python virtual environment:

    ```sh
    source .venv/bin/activate
    ```

2. Open the Jupyter Notebook server:

    ```sh
    jupyter notebook
    ```

3. Navigate to the notebook files (e.g., `ncbi-cell/ipynb/Chapter-01-CELLxGENE.ipynb`) and start exploring! The notebooks are intended to be processed in order, from chapter 01 through chapter 05.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
