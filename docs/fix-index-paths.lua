#!/usr/bin/env lua

-- fix-index-paths.lua
-- Post-render script to fix paths in the generated site
-- This script:
-- 1. Copies docs/index.html to the site root and fixes relative paths
-- 2. Copies markdown files from _site/docs back to docs directory

-- Helper function to check if file exists
local function file_exists(path)
	local f = io.open(path, "r")
	if f ~= nil then
		io.close(f)
		return true
	else
		return false
	end
end

-- Helper function to check if directory exists
local function dir_exists(path)
	-- Try to open the directory as a file, which will fail
	-- Then try using a platform-agnostic approach
	local f = io.open(path .. "/.", "r")
	if f then
		io.close(f)
		return true
	end
	return false
end

-- Helper function to read file
local function read_file(path)
	local f = io.open(path, "r")
	if not f then
		return nil
	end
	local content = f:read("*all")
	f:close()
	return content
end

-- Helper function to write file
local function write_file(path, content)
	local f = io.open(path, "w")
	if not f then
		return false
	end
	f:write(content)
	f:close()
	return true
end

-- Helper function to copy file
local function copy_file(src, dest)
	local content = read_file(src)
	if content then
		return write_file(dest, content)
	end
	return false
end

-- Helper function to list files in directory
local function list_files(dir, pattern)
	local files = {}
	local command

	-- Detect OS and use appropriate command
	if package.config:sub(1, 1) == "\\" then
		-- Windows
		command = 'dir /b "' .. dir .. '" 2>nul'
	else
		-- Unix-like (Linux, macOS)
		command = 'ls "' .. dir .. '" 2>/dev/null'
	end

	local handle = io.popen(command)
	if handle then
		for file in handle:lines() do
			if pattern and file:match(pattern) then
				table.insert(files, file)
			elseif not pattern then
				table.insert(files, file)
			end
		end
		handle:close()
	end

	return files
end

-- Main function
local function main()
	-- Define paths relative to project root
	local site_dir = "_site"
	local docs_dir = "docs"
	local source_index = site_dir .. "/docs/index.html"
	local dest_index = site_dir .. "/index.html"

	-- Part 1: Copy and fix index.html
	if file_exists(source_index) then
		local content = read_file(source_index)
		if content then
			-- Fix relative paths
			content = content:gsub('href="../docs/', 'href="docs/')
			content = content:gsub('src="../site_libs/', 'src="site_libs/')

			if write_file(dest_index, content) then
				print("✓ Fixed paths in root index.html")
			else
				print("⚠️  Failed to write fixed index.html")
			end
		else
			print("⚠️  Failed to read " .. source_index)
		end
	else
		print("⚠️  No _site directory found - using default homepage")
	end

	-- Part 2: Copy markdown files back to docs directory
	local site_docs_dir = site_dir .. "/docs"
	if dir_exists(site_docs_dir) then
		local md_files = list_files(site_docs_dir, "%.md$")
		local copied = false

		for _, file in ipairs(md_files) do
			local src = site_docs_dir .. "/" .. file
			local dest = docs_dir .. "/" .. file
			if copy_file(src, dest) then
				copied = true
			end
		end

		if copied then
			print("✓ Markdown files copied to docs/")
		else
			print("⚠️  No markdown files found")
		end
	end
end

-- Run the script
main()
