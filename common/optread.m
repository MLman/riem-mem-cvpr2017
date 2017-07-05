function v = optread(option, name, default)
if isfield(option, name)
    eval(sprintf('v=option.%s;',name));
else
    v = default;
end