
# coding: utf-8

# In[324]:


import numpy as np

#Creamos los arrays con los valores de Ltot y Rtot sobre los cuales vamos a obtener su error

lums = np.arange(22.375, 22.396, 0.001)
rads = np.arange(10.975, 10.99, 0.001)
error_mat = []

#Evaluamos el codigo original

for Rtot in rads:
    for Ltot in lums:
        
        get_ipython().magic('run Codigo_auxiliar.ipynb')
        
        error_mat.append(Err_rel_total)
        
        print(Ltot, Rtot, Tcentral, Err_rel_total)

error_mat = np.array(error_mat).reshape(len(rads), len(lums))
error_mat


# In[325]:


lums = np.arange(22.375, 22.396, 0.001)
rads = np.arange(10.975, 10.99, 0.001)


# In[326]:


lums = lums.round(1)
rads = rads.round(1)


# In[327]:


lums = np.array(lums).tolist()
rads = np.array(rads).tolist()
type(error_mat)


# In[328]:


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Display matrix

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }

sns.set()

plt.figure(figsize=(40,21))
resultado = sns.heatmap(error_mat, cmap = 'hot', linewidths= 0.01, square=True)
#resultado = sns.heatmap(error_mat, cmap = 'coolwarm', annot = error_mat, linewidths= 0.5)

resultado.set_xticks(np.arange(len(lums)))
resultado.set_yticks(np.arange(len(rads)))

resultado.set_xticklabels(lums, fontdict=font)
resultado.set_yticklabels(rads, fontdict=font)

plt.title('Estudio del error relativo total', fontsize=40 )
plt.xlabel('luminosidad total ($L_{tot}$)', fontsize=35)
plt.ylabel('radio total ($R_{tot}$)', fontsize=35)

plt.savefig('mapacalorgrande.png')
plt.show(resultado)


# In[329]:


error_minimo = error_mat.min()
error_minimo


# In[330]:


np.save('error_mat', error_mat)


# ### Plot formatting ...

# In[331]:


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')


# In[332]:


lums = np.arange(22.375, 22.396, 0.001).round(3)
rads = np.arange(10.975, 10.99, 0.001).round(3)
print(lums)
print(rads)


# In[333]:


import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
plt.rcParams["font.family"] = "Times New Roman"


# In[334]:


im = plt.imshow(error_mat, cmap='nipy_spectral');


# In[335]:


import seaborn as sns

fig = plt.figure(figsize=(20,8))
ax = sns.heatmap(error_mat, cmap = 'nipy_spectral', linewidths= 0.01, cbar=False, square=True);

ax.set_xlabel('\nluminosidad total (L$_{tot}$)', fontsize=24);
ax.set_ylabel('radio total (R$_{tot}$)\n', fontsize=24);

lum_label = [];
for i in np.arange(0, len(lums), 1):
    if i%2 == 0:
        lum_label.append(lums[i])
    else:
        lum_label.append(' ');
        
rad_label = [];
for i in np.arange(0, len(rads), 1):
    if i%2 == 0:
        rad_label.append(rads[i])
    else:
        rad_label.append(' ');
    
ax.set_xticklabels(lum_label);
ax.set_yticklabels(rad_label);
ax.tick_params(labelsize = 18, axis = 'x', which = 'major', pad=10, direction='inout', length = 6, width = 1);
ax.tick_params(labelsize = 18, axis = 'y', which = 'major', pad=10, direction='inout', length = 6, width = 1);
    
plt.title('Estudio del error relativo total\n', fontsize=32);
plt.xticks(fontsize=18, rotation=0);
plt.yticks(fontsize=18, rotation=0);

from mpl_toolkits.axes_grid1 import make_axes_locatable
color_code = np.array([]);
divider = make_axes_locatable(ax);
cax = divider.append_axes('right', size='3%', pad=0.1);
cb = fig.colorbar(im, cax=cax);
#cb.set_label('relative error', fontsize=16, rotation=-90);
cb.solids.set_edgecolor('face');
cax.tick_params(labelsize=14, axis='y', which='major', pad=10, direction='inout', length=4, width=1);

#ax.annotate('Error mínimo', fontsize = 20,
            #xy=(4,13), xycoords='data',
            #xytext=(-180, -260), textcoords='offset points',
            #arrowprops=dict(arrowstyle="->"));


fig.tight_layout();
fig.savefig('rel_error_heatmap.png', bbox_inches='tight');
plt.show();


# In[342]:


print(error_minimo)
print(lums[7])
print(rads[7])
error_mat[7,7]

