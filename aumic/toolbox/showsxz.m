function showsxz(model)

nz = size(model.z,2);
ifam = 1;
zname{1,1} = 'EM';
for iz = 1:nz
    if iz > sum(model.nzInFam(1:ifam))
        ifam = ifam+1;
    end
    zname{1,iz+1} = ['z',num2str(iz),'/f',num2str(ifam)];
end
xsxz = [zname;model.x num2cell(model.sxz)];
disp(xsxz)