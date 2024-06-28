
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

zeropoint_dir = {
    'F115W': {
        'CEERS-NE': 28.9,
        'CEERS-SW': 28.9,
    },
              
    'F160W': {
        'CEERS-NE': 28.9,
        'CEERS-SW': 28.9,
        'AEGIS': 25.946,
    },
}

run_strs = [
    'F115W_single',
    'F115W_multi',
    'F160W_single',
    'F160W_multi',
    'F160W_vdw',
    ]

fi_dir = {
    'F115W_single': 'F115W',
    'F115W_multi': 'F115W',
    'F160W_single': 'F160W',
    'F160W_multi': 'F160W',
    'F160W_vdw': 'F160W',
    }

how_dir = {
    'F115W_single': 0,
    'F160W_single': 0,
    #'F115W_multi': 0,
    #'F160W_multi': 0,
    'F160W_vdw': 0,
    }

how_param_dir = {
    'F115W_single': [1, 1, 1, 1, 1, 1, 1],
    'F160W_single': [1, 1, 1, 1, 1, 1, 1],
    #'F115W_multi': [1, 1, 1, 1, 1, 1, 1],
    #'F160W_multi': [1, 1, 1, 1, 1, 1, 1],
    'F160W_vdw': [1, 1, 0, 0, 0, 0, 0],
}

ix_tofit = [
    #1509,
    #3788,
    #4000,
    #4768,
    #719,
    #1743,
    #4317,
    #3928,
    ]

constraints_dir = {
    'F115W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5',
    'F160W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3',
    'F115W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F160W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F160W_vdw': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3',
    }
