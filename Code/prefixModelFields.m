function model = prefixModelFields(model, prefix)
    model.mets = strcat(prefix, model.mets);
    model.rxns = strcat(prefix, model.rxns);
    model.genes = strcat(prefix, model.genes);
end