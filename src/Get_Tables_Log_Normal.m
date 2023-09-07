% cd /Users/cangtao/cloud/GitHub/NANOGrav-SIGW/src
cd /Users/cangtao/cloud/Matlab/DM_BH/SIGW
clear

reload = 1

s = logspace(-2, 1, 200)
x = logspace(-5, 5, 100)

s = logspace(-2, 1, 400)
x = logspace(-5, 5, 200)

H5_File = '/Users/cangtao/cloud/GitHub/NANOGrav-SIGW/data/GW_Tables.h5'

clc

ns = length(s);
nx = length(x);

Estimated_Time = ns*nx*0.1/60

if reload
    tic
    for sid = 1:ns
        for xid = 1:nx
            s_ = s(sid);
            x_ = x(xid);

            r = OmegaGW(1, s_, 1, x_, 2);
            Tab(sid, xid) = r;

        end
        status = sid/ns
        toc
    end

    save data/NanoGrav_Tab.mat Tab
else
    load data/NanoGrav_Tab.mat
end


delete(H5_File)
TabSize=size(Tab);
h5create(H5_File, '/Sigma_axis', ns);
h5create(H5_File, '/x_axis', nx);
h5create(H5_File, '/dGWdlnk_Tables', TabSize);

h5write(H5_File, '/Sigma_axis', s);
h5write(H5_File, '/x_axis', x);
h5write(H5_File, '/dGWdlnk_Tables', Tab);
