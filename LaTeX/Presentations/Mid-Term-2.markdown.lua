local ran_ok, error = pcall(function() local kpse = require("kpse") kpse.set_program_name("luatex") local lfs = require("lfs") local cacheDir = "./_markdown_Mid-Term-2" if not lfs.isdir(cacheDir) then assert(lfs.mkdir(cacheDir)) end local md = require("markdown") local convert = md.new({cacheDir = "./_markdown_Mid-Term-2", citations = true, definitionLists = true, footnotes = true, hashEnumerators = true, hybrid = true, pipeTables = true, smartEllipses = true, tableCaptions = true, tightLists = false, } ) local input = assert(io.open("./Mid-Term-2.markdown.in", "r"):read("*a")) print(convert(input:gsub("\r\n?", "\n"))) end) if not ran_ok then local file = io.open("./Mid-Term-2.markdown.err", "w") if file then file:write(error .. "\n") file:close() end print('\\markdownError{An error was encountered while executing Lua code}{For further clues, examine the file "./Mid-Term-2.markdown.err"}') end
