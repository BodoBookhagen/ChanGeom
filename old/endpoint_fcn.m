function is_end_point = endpoint_fcn(nhood)
% function is_end_point = endpoint_fcn(nhood)

interval1 = [0 1 0; -1 1 -1; -1 -1 -1];
interval2 = [1 -1 -1; -1 1 -1; -1 -1 -1];

for k= 1:4
    C = bwhitmiss(nhood, rot90(interval1,k));
    D = bwhitmiss(nhood, rot90(interval2,k));
    if (C(2,2) == 1) || (D(2,2) == 1)
        is_end_point = true;
        return
    end
end
is_end_point = false;