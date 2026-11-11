-- Usage: osascript render_docx.applescript input.docx output.pdf
on run argv
    set inputFile to POSIX file (item 1 of argv)
    set outputFile to POSIX file (item 2 of argv)
    tell application "Pages"
        set theDocument to open inputFile
        delay 3
        export theDocument to outputFile as PDF
        close theDocument saving no
    end tell
    return item 2 of argv
end run
