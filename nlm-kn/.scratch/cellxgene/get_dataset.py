#!/usr/bin/env python
# coding: utf-8

# # Fetch full metadata for a Dataset
# See: https://github.com/chanzuckerberg/single-cell-curation/blob/0c77179d2e794846861f8109c037b723507959cb/notebooks/curation_api/python_raw/get_dataset.ipynb

# The script in this notebook retrieves full metadata for a given Dataset.
# 
# Fetching a Dataset requires only the Collection id and the Dataset id; it does not require an API key/access token.

# ### Import dependencies

# In[1]:


import requests


# #### <font color='#bc00b0'>Please fill in the required values:</font>

# <font color='#bc00b0'>(Required) Enter the id of the Collection that contains the Dataset for which you want to fetch full metadata</font>
# 
# _The Collection id can be found by looking at the url path in the address bar 
# when viewing your Collection in the CZ CELLxGENE Discover data portal: `/collections/{collection_id}`._

# In[7]:


collection_id = "28e9d721-6816-48a2-8d0b-43bf0b0c0ebc"


# <font color='#bc00b0'>(Required) Enter the id of the Dataset</font>
# 
# _The Dataset id can be found by using the `/collections/{collection_id}` endpoint and filtering for the Dataset of interest OR by looking at the url path in the address when viewing your Dataset using the CZ CELLxGENE Explorer browser tool: `/e/{dataset_id}.cxg/`._

# In[8]:


dataset_id = "6e00ccf7-0749-46ef-a999-dba785630d52"


# ### Specify domain (and API url)

# In[9]:


domain_name = "cellxgene.cziscience.com"
site_url = f"https://{domain_name}"
api_url_base = f"https://api.{domain_name}"


# ### Formulate request and fetch a Datasets metadata

# In[10]:


dataset_path = f"/curation/v1/collections/{collection_id}/datasets/{dataset_id}"
url = f"{api_url_base}{dataset_path}"
res = requests.get(url=url)
res.raise_for_status()
res_content = res.json()
print(res_content)


# ### Download Dataset Assets
# 
# The dataset metadata provides download URLs for every asset associated with the current dataset version. For public collections, that means the most recently published version of a dataset. For private collections, that means the most recently successfully processed dataset version.
# 
# These download URLs are permalinks to download the assets for this particular version of a dataset. If this dataset is revised, you would need to fetch the dataset metadata again to get the latest dataset version asset download links.

# In[11]:


assets = res_content["assets"]
dataset_id = res_content["dataset_id"]
for asset in assets:
    download_filename = f"{dataset_id}.{asset['filetype']}"
    print(f"\nDownloading {download_filename}... ")
    with requests.get(asset["url"], stream=True) as res:
        res.raise_for_status()
        filesize = int(res.headers["Content-Length"])
        with open(download_filename, "wb") as df:
            total_bytes_received = 0
            for chunk in res.iter_content(chunk_size=1024 * 1024):
                df.write(chunk)
                total_bytes_received += len(chunk)
                percent_of_total_upload = float("{:.1f}".format(total_bytes_received / filesize * 100))
                color = "\033[38;5;10m" if percent_of_total_upload == 100 else ""
                print(f"\033[1m{color}{percent_of_total_upload}% downloaded\033[0m\r", end="")
print("\n\nDone downloading assets")


# In[ ]:




