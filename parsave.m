function parsave(filename, vname, v)
eval([vname '= v;'])
save(filename, vname, '-append')
end