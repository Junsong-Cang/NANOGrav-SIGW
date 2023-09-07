cd /Users/cangtao/cloud/Matlab/DM_BH/SIGW
clear

reload = 0
H5_File = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/GW_Tables.h5'
A = 1e8

freq_range = [-9, -7]
freq_bh_range = [-9.5, -6.5]
Sigma_range = [-2, 1]

nk = 200
nkbh = 200
ns = 100

nk = 200
nkbh = 250
ns = 300

clc
const = log10(1.546e-15);
k1 = freq_range(1)-const;
k2 = freq_range(2)-const;

kbh1 = freq_bh_range(1)-const
kbh2 = freq_bh_range(2)-const

! cp ~/FUNCTIONS/linspace.m ./
ks= linspace(k1, k2, nk);
kbhs= linspace(kbh1, kbh2, nkbh);
Sigmas = linspace(Sigma_range(1), Sigma_range(2), ns);
! rm linspace.m

ks = 10.^ks;
kbhs = 10.^kbhs;
Sigmas = 10.^Sigmas;


if reload
    tic
    DoF = Find_DoF(ks);
    for sid = 1:ns
        for kbhid = 1:nkbh

            [k_axis, GW_axis] = OmGW_Interp(A, Sigmas(sid), kbhs(kbhid));

            for kid = 1:nk

                DoF_here = DoF(kid);
                
                dOmGW_dlnk = interp1(k_axis, GW_axis, ks(kid));
                dOmGW_dlnk0 = 0.38*9.1E-5*((DoF_here/106.75)^(-1/3)).*dOmGW_dlnk;
                Tab(sid, kbhid, kid) = dOmGW_dlnk0;

            end

        end
        status = sid/ns
    end
    toc
    save data/NanoGrav_Tab.mat Tab
else
    load data/NanoGrav_Tab.mat
end

freq_axis = 1.546e-15*ks;

delete(H5_File)
TabSize=size(Tab);
h5create(H5_File, '/Sigma_axis', ns);
h5create(H5_File, '/freq_axis', nk);
h5create(H5_File, '/kbh_axis', nkbh);
h5create(H5_File, '/dGWdlnk_Tables', TabSize);

h5write(H5_File, '/Sigma_axis', Sigmas);
h5write(H5_File, '/freq_axis', freq_axis);
h5write(H5_File, '/kbh_axis', kbhs);
h5write(H5_File, '/dGWdlnk_Tables', Tab);
