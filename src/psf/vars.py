

data_path = '/Users/rschmidt-eder/cassi2024/data'
pdata_path = '/Users/rschmidt-eder/cassi2024/processed_data/psf'
res_path = '/Users/rschmidt-eder/cassi2024/results/jwst_goods-s_size'


filters = [
    'F115W',
    #'F160W',
]

fields = ['GOODS-S',
          #'CEERS-SW',
          #'AEGIS',
          ]


sci_dir = {
    'F115W': {
        'GOODS-S': '{}/GOODS-S_mosaic/gds-grizli-v7.2-f115w-clear_drc_sci.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f115w-clear_drc_sci.fits'.format(data_path),
    },
              
    'F160W': {
        'CEERS-NE': '{}/CEERS-NE/ceers-ne-grizli-v4.0-f160w_drz_sci.fits'.format(data_path),
        'CEERS-SW': '{}/CEERS-SW/ceers-sw-grizli-v4.0-f160w_drz_sci.fits'.format(data_path),
        'AEGIS': '{}/AEGIS/gbrammer/egs-100mas-f160w_drz_sci.fits'.format(data_path),
    },
}


starlist_dir = {
    'F115W': {
        'GOODS-S': '{}/starlists/starlist_GOODS-S.xy'.format(pdata_path),
        'CEERS-SW': '{}/starlists/starlist_CEERS-SW.xy'.format(pdata_path),
    },
    'F160W': {
        'CEERS-NE': '{}/starlists/starlist_F160W_CEERS-NE.xy'.format(pdata_path),
        'CEERS-SW': '{}/starlists/starlist_F160W_CEERS-SW.xy'.format(pdata_path),
        'AEGIS': '{}/starlists/starlist_F160W_AEGIS.xy'.format(pdata_path)
    },
}

psfin_dir = {
    'F115W': '{}/GOODS-S_mosaic/psf-good-s-nrca3-f115w.fits'.format(data_path), #0.04"/px
    'F160W': '{}/CEERS-NE/psf-ceers-f160w.fits'.format(data_path) #0.04"/px
    #'F160W': '{}/AEGIS/aegis_3dhst_v4.0_wfc3_psf/aegis_3dhst.v4.0.F160W_psf.fits'.format(data_path), #from 3DHST. better use Gabe's if I'm using Gabe's sci mosaic. 0.06"/px
}

pxscalein_dir = {
    'F115W': 0.04,
    'F160W': 0.04,
}

pxscale_mosaic_dir = {
    'F115W': [0.02],
    'F160W': [0.04,
              0.1],
}
    
