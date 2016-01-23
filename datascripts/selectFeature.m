function feature_table = selectFeature( data, feature_name )
    feature_indicator = cellfun(@(x)( ~isempty(x) ), regexp(data.ID, feature_name));
    feature_table = data(feature_indicator,:);
end