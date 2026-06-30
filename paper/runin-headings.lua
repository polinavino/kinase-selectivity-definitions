-- achemso environments that pandoc keeps as untitled Divs; add their section
-- heading so they appear as titled sections in the docx, matching the PDF.
local env_titles = {
  acknowledgement = "Acknowledgement",
  suppinfo = "Supporting Information Available",
}

function Div(div)
  for cls, title in pairs(env_titles) do
    if div.classes:includes(cls) then
      local blocks = { pandoc.Header(1, pandoc.Str(title)) }
      for _, b in ipairs(div.content) do
        table.insert(blocks, b)
      end
      return blocks
    end
  end
  return nil
end

-- Number figures and tables so captions read "Figure N." / "Table N." in
-- document order, matching the PDF (pandoc does not auto-number for docx).
local fig_count = 0
local tbl_count = 0

local function prepend_label(caption, label)
  if caption and caption.long and #caption.long > 0 then
    local first = caption.long[1]
    if first.content then
      table.insert(first.content, 1, pandoc.Space())
      table.insert(first.content, 1, pandoc.Str(label))
    end
  end
end

function Figure(fig)
  fig_count = fig_count + 1
  prepend_label(fig.caption, "Figure " .. fig_count .. ".")
  return fig
end

function Table(tbl)
  tbl_count = tbl_count + 1
  prepend_label(tbl.caption, "Table " .. tbl_count .. ".")
  return tbl
end

-- Convert LaTeX \paragraph run-in headings (pandoc Header level >= 4) into
-- run-in bold text merged into the following paragraph, matching the PDF where
-- the heading and its body start on the same line.
function Pandoc(doc)
  local blocks = doc.blocks
  local out = {}
  local i = 1
  while i <= #blocks do
    local b = blocks[i]
    if b.t == "Header" and b.level >= 4 then
      local title = pandoc.Strong(b.content)
      local nxt = blocks[i + 1]
      if nxt and nxt.t == "Para" then
        local inlines = { title, pandoc.Space() }
        for _, inl in ipairs(nxt.content) do
          table.insert(inlines, inl)
        end
        table.insert(out, pandoc.Para(inlines))
        i = i + 2
      else
        table.insert(out, pandoc.Para({ title }))
        i = i + 1
      end
    else
      table.insert(out, b)
      i = i + 1
    end
  end
  return pandoc.Pandoc(out, doc.meta)
end
