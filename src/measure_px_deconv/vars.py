# I belive this file should be the same as that in /src/galfit/vars.py
# So could copy that file into /src/size/ and replace this one



import os
home = os.path.expanduser('~')

data_path = '{}/cassi2024/data'.format(home)
pdata_path = '{}/cassi2024/processed_data'.format(home)
res_path = '{}/cassi2024/results'.format(home)

filters = [
    'F115W',
    'F277W',
    'F356W',
    'F444W',
]

field_ord = [
    'GOODS-S',
    ]

pxsc_dir = {
    'F115W': {
        'GOODS-S': 0.03121,
    },
    
    'F277W': {
        'GOODS-S': 0.06294,
    },
    
    'F356W': {
        'GOODS-S': 0.06294,
    },
    
    'F444W': {
        'GOODS-S': 0.06294,
    },
}

zeropoint_dir = {
    'F115W': {
        'GOODS-S': 28.022794043870377,
    },
    
    'F277W': {
        'GOODS-S': 26.489416915579795,
    },
    
    'F356W': {
        'GOODS-S': 26.489416915579795,
    },
    
    'F444W': {
        'GOODS-S': 26.489416915579795,
    },
    #
    #'F115W': {
    #    'CEERS-NE': 28.9,
    #    'CEERS-SW': 28.9,
    #},
    #          
    #'F160W': {
    #    'CEERS-NE': 28.9,
    #    'CEERS-SW': 28.9,
    #    'AEGIS': 25.946,
    #},
}

run_strs = [
    'F115W_single',
    'F115W_multi',
    'F277W_single',
    'F277W_multi',
    'F356W_single',
    'F356W_multi',
    'F444W_single',
    'F444W_multi',
    #'F160W_single',
    #'F160W_multi',
    #'F160W_vdw',
    ]

fi_dir = {
    'F115W_single': 'F115W',
    'F115W_multi': 'F115W',
    'F277W_single': 'F277W',
    'F277W_multi': 'F277W',
    'F356W_single': 'F356W',
    'F356W_multi': 'F356W',
    'F444W_single': 'F444W',
    'F444W_multi': 'F444W',
    'F160W_single': 'F160W',
    'F160W_multi': 'F160W',
    'F160W_vdw': 'F160W',
    }

how_dir = {
    'F115W_single': 0,
    'F277W_single': 0,
    'F356W_single': 0,
    'F444W_single': 0,
    'F160W_single': 0,
    #'F115W_multi': 0,
    #'F160W_multi': 0,
    'F160W_vdw': 0,
    }

how_param_dir = {
    'F115W_single': [1, 1, 1, 1, 1, 1, 1],
    'F277W_single': [1, 1, 1, 1, 1, 1, 1],
    'F356W_single': [1, 1, 1, 1, 1, 1, 1],
    'F444W_single': [1, 1, 1, 1, 1, 1, 1],
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
    'F277W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5',
    'F356W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5',
    'F444W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5',
    'F160W_single': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3',
    'F115W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F277W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F356W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F444W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 5\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F160W_multi': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3\n2 x -10 10\n2 y -10 10\n2 n 0.2 to 8\n2 q 0.0001 to 1\n2 mag -3 5',
    'F160W_vdw': '1 x -10 10\n1 y -10 10\n1 n 0.2 to 8\n1 q 0.0001 to 1\n1 mag -3 3',
    }
