# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.4
#   kernelspec:
#     display_name: Python [conda env:bbsnippets]
#     language: python
#     name: conda-env-bbsnippets-py
# ---

# %%
import h5py
import re
import pandas as pd
import numpy as np
import pathlib


# %% [markdown]
# # Examples how to extract information from ilastik projects
#
# This shows how to:
# - Export a csv that can be loaded as metata csv in cellprofiler to reproduce the 'random' cropping
# - How to dump the labels used for training
#
#
# Further useful scripts how to use ilastik from the command line, e.g. how to train it with pre-define labels, can be found here: https://github.com/ilastik/ilastik/tree/master/bin
#
# By
# vito.zanotelli@gmail.com

# %%
class config:
    fn_ilastik = './exampledata/example_training.ilp'
    fn_crop_list ='~/tmp/crop_metadata.csv' # This can be used within cellprofiler to reproduce the croping identically to the training data
    fol_out_labels = pathlib.Path('/home/vitoz/tmp')
C=config

class variables:
    re_crop = re.compile('(?P<basename>.*)(_s(P<scale>[0-9]+))?_x(?P<x>[0-9]+)_y(?P<y>[0-9]+)_w(?P<w>[0-9]+)_h(?P<h>[0-9]+).*')
    COL_FN = 'filename_training'
    SUFFIX_LABEL = '_label.npy'
V=variables


# %% [markdown]
# Helper function to dump a `hdf5` content nicely formated from:
# https://stackoverflow.com/a/53340677

# %%
def descend_obj(obj,sep='\t'):
    """
    Iterate through groups in a HDF5 file and prints the groups and datasets names and datasets attributes
    """
    if type(obj) in [h5py._hl.group.Group,h5py._hl.files.File]:
        for key in obj.keys():
            print(sep,'-',key,':',obj[key])
            descend_obj(obj[key],sep=sep+'\t')
    elif type(obj)==h5py._hl.dataset.Dataset:
        for key in obj.attrs.keys():
            print(sep+'\t','-',key,':',obj.attrs[key])

def h5dump(path,group='/'):
    """
    print HDF5 file metadata

    group: you can give a specific group, defaults to the root group
    """
    with h5py.File(path,'r') as f:
         descend_obj(f[group])


# %% [raw]
# h5dump(C.fn_ilastik)

# %% [markdown]
# ## Dump a metadata csv to reproduce the random crops
#
# Make a metadata file from the ilasik project filenames, that can be used to reproduce the crops by using this csv as metadta file.
# This file can be loaded in cellprofiler as metadata. The `crop bb` module from ImcPluginsCP (https://github.com/BodenmillerGroup/ImcPluginsCP) can use metadata as parameters to specify were to crop.

# %%
with h5py.File(C.fn_ilastik, 'r') as f:
    lanes = f['Input Data']['infos'].values()
    fns_training = [lane['Raw Data/nickname'][()].decode('UTF-8') for lane in lanes ]


# %%
def split_names(x, re_comp):
    c = re_comp
    m = c.match(x)
    g = m.groups() 
    return pd.Series({l: g[i-1] for l, i in c.groupindex.items()}, name=x.index)


# %%
dat_training = pd.DataFrame({V.COL_FN: fns_training})
dat_training = dat_training.join(dat_training[V.COL_FN].apply(split_names, re_comp=V.re_crop))

# %%
dat_training

# %%
dat_training.to_csv(C.fn_crop_list, index=False)

# %% [markdown]
# ## Extract ilastik training labels
#
# These labels will be saved in the label output folder as a `.npy` array.
#
# This is usefull e.g. to combine classifiers.
#

# %%
with h5py.File(C.fn_ilastik, 'r') as f:
    labels = f['/PixelClassification/LabelSets']
    lanes = f['Input Data']['infos']
    for label, lane in zip(labels.values(),lanes.values()):
        name = lane['Raw Data/nickname'][()].decode('UTF-8')
        for val in label.values():
            print(name)
            np.save(C.fol_out_labels / (name+V.SUFFIX_LABEL), val[:])

# %% [markdown]
# # Environment

# %%
import sys
# !conda env export -p {sys.prefix}
