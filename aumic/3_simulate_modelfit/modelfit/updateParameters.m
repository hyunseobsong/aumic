function kinetic = updateParameters(kinetic,para,paraopt)

kmax.num = cell2mat(kinetic.kmax(:,2));
kmax.name = kinetic.kmax(:,1);
K.num = cell2mat(kinetic.K(:,2));
K.name = kinetic.K(:,1);

% search kmax
np = length(para);
iip.kmax = []; ip_ = [];
for ip = 1:np
    ip_ = find(strcmp(para{ip},kmax.name));
    if ~isempty(ip_)
        kinetic.kmax{ip_,2} = paraopt(ip);
%         kinetic.kmax{ip_,2} = char(num2str(paraopt(ip)));
%         kinetic.kmax(ip_,2) = num2str(paraopt(ip));
    end
%     iip.kmax = [iip.kmax;ip_];
end

% search K
iip.K= []; ip_ = [];
for ip = 1:np
    ip_ = find(strcmp(para{ip},K.name));
    if ~isempty(ip_)
        kinetic.K{ip_,2} = paraopt(ip);
%         kinetic.K{ip_,2} = char(num2str(paraopt(ip)));
    end
%     iip.K = [iip.K;ip_];
end
