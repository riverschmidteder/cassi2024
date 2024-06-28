

data_path = '/home/tonypahl/jwst_ceers_size/data'
pdata_path = '/home/tonypahl/jwst_ceers_size/processed_data'
res_path = '/home/tonypahl/jwst_ceers_size/results'


filters = [
    'F115W',
    'F160W',
]

field_ord = [
    'CEERS-NE',
    'CEERS-SW',
    'AEGIS',
    ]

sci_dir = {
    'F115W': {
        'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-f115w-clear_drc_sci.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f115w-clear_drc_sci.fits'.format(data_path),
    },
              
    'F160W': {
        'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-f160w_drz_sci.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f160w_drz_sci.fits'.format(data_path),
        'AEGIS': '{}/AEGIS/gbrammer/egs-100mas-f160w_drz_sci.fits'.format(data_path),
    },
}

wht_dir = {
    'F115W': {
        'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-f115w-clear_drc_wht.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f115w-clear_drc_wht.fits'.format(data_path),
    },
              
    'F160W': {
        'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-f160w_drz_wht.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f160w_drz_wht.fits'.format(data_path),
        'AEGIS': '{}/AEGIS/gbrammer/egs-100mas-f160w_drz_wht.fits'.format(data_path),
    },
}

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


pxsc_dir = {
    'F115W': {
        'CEERS-NE': 0.02,
        'CEERS-SW': 0.02,
    },
              
    'F160W': {
        'CEERS-NE': 0.04,
        'CEERS-SW': 0.04,
        'AEGIS': 0.1,
    },
}
