function text = describeFileContent(fileContent, options)
%DESCRIBEFILECONTENT Format an importAgilentUV/CH/MS output struct as multiline text.
%
% Syntax:
%   text = describeFileContent(fileContent)
%   text = describeFileContent(fileContent, ExcludeFields="method.datetime")
%
% Description:
%   Turns the metadata returned by importAgilentUV, importAgilentCH or importAgilentMS into a
%   human-readable string, e.g. for logging or for a free-text description field that has no
%   dedicated property of its own.
%   Every top-level substruct (file, sample, method, instrument, signal, ...) is included by
%   default, one line per field, formatted as "structName.fieldName: value". Fields that are not
%   scalar or cannot be converted to a single string (e.g. fileContent.signal, which holds
%   vectors/matrices) are skipped automatically.
%
% Input Arguments:
%   fileContent - (1,1) struct
%       A struct as returned by importAgilentUV, importAgilentCH or importAgilentMS (or any
%       struct with the same shape: top-level substructs of scalar fields).
%
% Name-Value Arguments:
%   ExcludeFields - (1,:) string
%       Additional fields to leave out, on top of the non-scalar fields skipped automatically.
%       Write "structName.fieldName" to exclude a single field, or just "structName" to exclude an
%       entire top-level substruct. Naming a field or substruct that is not present, or that is
%       already skipped automatically, is not an error. Default: string.empty (nothing
%       additionally excluded).
%
% Output Arguments:
%   text - (1,1) string
%       The formatted description, one field per line, joined with newline.

arguments
    fileContent (1,1) struct
    options.ExcludeFields (1,:) string = string.empty
end

excludedStructNames = options.ExcludeFields(~contains(options.ExcludeFields, "."));
excludedFieldPaths = options.ExcludeFields(contains(options.ExcludeFields, "."));

structNames = setdiff(string(fieldnames(fileContent))', excludedStructNames, "stable");

lines = string.empty(0, 1);
for structName = structNames
    lines = [lines; localStructLines(structName, fileContent.(structName), excludedFieldPaths)]; %#ok<AGROW>
end

text = join(lines, newline);

end

function lines = localStructLines(structName, s, excludedFieldPaths)
% Format every field of struct s as "structName.fieldName: value" lines, skipping the fields named
% in excludedFieldPaths and any field that is not scalar or cannot be converted to a single string
% (e.g. a vector/matrix field).

excludedFieldNames = extractAfter(excludedFieldPaths(startsWith(excludedFieldPaths, structName + ".")), structName + ".");
fieldNames = setdiff(string(fieldnames(s))', excludedFieldNames, "stable");

lines = string.empty(0, 1);
for name = fieldNames
    value = s.(name);
    if ~isscalar(value)
        continue
    end
    try
        valueText = string(value);
    catch
        continue
    end
    if ~isscalar(valueText)
        continue
    end
    lines(end+1, 1) = structName + "." + name + ": " + valueText; %#ok<AGROW>
end

end
