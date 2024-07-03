
#Path to the data folder (not a specific file/folder)
import os
home = os.path.expanduser('~')
data_path = '{}/cassi2024/data'.format(home)
pdata_path = '{}/cassi2024/processed_data'.format(home)
res_path = '{}/cassi2024/results'.format(home)

#The filters we're cutting stamps out of
filters = [
    'F115W',
    #'F160W',
]
#Not that important
field_ord = [
    #'CEERS-NE',
    #'CEERS-SW',
    #'AEGIS',
    'GOODS-S',
    'GOODS-N',
    ]

#This is for the science image
sci_dir = {
    'F115W': {
        'GOODS-N': '{}/GOODS-N_mosaic/gdn-grizli-v7.3-f115w-clear_drc_sci.fits'.format(data_path),
        'GOODS-S': '{}/GOODS-S_mosaic/gds-grizli-v7.2-f115w-clear_drc_sci.fits'.format(data_path),
    },

}

#This is for the associated weight image
wht_dir = {
    'F115W': {
        'GOODS-N': '{}/GOODS-N_mosaic/gdn-grizli-v7.3-f115w-clear_drc_wht.fits'.format(data_path),
        'GOODS-S': '{}/GOODS-S_mosaic/gds-grizli-v7.2-f115w-clear_drc_wht.fits'.format(data_path),
    },
}

#Once the segmentation image is made we'll put it here. That's what the source extractor makes.
seg_dir = {
    'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-ir_seg.fits'.format(data_path),
    'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-ir_seg.fits'.format(data_path),
    'AEGIS': '{}/AEGIS/gbrammer/egs-100mas-f160w_drz_seg.fits'.format(data_path),

    #'F115W': {
    #    'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-ir_seg.fits'.format(data_path),
    #    'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-ir_seg.fits'.format(data_path),
    #},
    #          
    #'F160W': {
    #    'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-ir_seg.fits'.format(data_path),
    #    'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-ir_seg.fits'.format(data_path),
    #    'AEGIS': '{}/AEGIS/gbrammer/egs-100mas-f160w_drz_seg.fits'.format(data_path),
    #},
}

#Defining the pixel scale of the mosaic
pxsc_dir = {
    'F115W': {
        'CEERS-NE': 0.02,
        'CEERS-SW': 0.02,
        'GOODS-N': 0.02,
        'GOODS-S': 0.02,
    },
              
    'F160W': {
        'CEERS-NE': 0.04,
        'CEERS-SW': 0.04,
        'AEGIS': 0.1,
    },
}
