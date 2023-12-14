#import "@preview/i-figured:0.2.2"
#import "@preview/physica:0.8.1": *
#import "@preview/pinit:0.1.2": *
#import "@preview/cetz:0.1.2"


#let conf(
  // 文章基础信息设置
  title: "Paper Title",
  subtitle: none,
  shortitle: none,
  author: none,
  email: none,
  date: datetime.today(),
  body

) = {
  // 字体信息设置，这里可以随意修改为自己想要使用的字体
  let song = "FZShuSong-Z01"
  let hei = "FZHei-B01"
  let kai = "FZKai-Z03"
  let xbsong = "FZXiaoBiaoSong-B05"
  let code = "DejaVu Sans Mono"

  let title-font = hei
  let author-font = kai
  let body-font = song
  let heading-font = xbsong
  let caption-font = kai
  let header-font = kai
  let strong-font = hei
  let emph-font = kai
  let raw-font = code

  // 设置文档元数据
  if {author != none} {
    set document(author: author, title: title)
  }

  // 在第一页的页眉添加标题
  set page(numbering: "1", number-align: center, header: align(left)[
    #set text(font: header-font)
    #shortitle
  ])

  // 设置英文、中文字体
  set text(font: ("Times New Roman", body-font))

  set heading(numbering: "1-1")
  set bibliography(style: "mla")
  
  show heading: it => box(width: 100%)[
    #v(0.50em)
    #set text(font: heading-font)
    #if it.numbering != none { counter(heading).display() }
    #h(0.75em)
    #it.body
]

  show heading.where(
    level: 1
  ): it => box(width: 100%)[
    #v(0.5em)
    #set align(left)
    #set heading(numbering: "一、")
    #it
    #v(0.75em)
  ]

  // Title
  align(center)[
    #block(text(font: title-font, weight: 700, 1.75em, title))
    #v(0.5em)
  ]

  // Author information.
  align(center)[
  #pad(
    top: 0.25em,
    bottom: 0.25em,
    x: 0.25em,
  )[
    #set text(font: ("Times New Roman", author-font))
    #author
    #par()[#email]
  ]
  ]

  // Main body
  set par(first-line-indent: 2em)
  set enum(indent: 2em)
  //set figure(gap: 0.8cm)
  //set math.equation(numbering: "(1)")

  show figure: it => [
    //#v(12pt)
    #set text(font: caption-font)
    #it
    #par()[#text(size: 0.0em)[#h(0.0em)]]
    //#v(12pt)
  ]

  show image: it => [
    #it
    #par()[#text(size: 0.0em)[#h(0.0em)]]
  ]

  show table: it => [
    #set text(font: body-font)
    #it
  ]

  show strong: set text(font: strong-font)
  show emph: set text(font: emph-font)
  show ref: set text(blue)
  show footnote: set text(blue)

  show raw.where(block: true): block.with(
    width: 100%,
    fill: luma(240),
    inset: 10pt,
  )

  show raw: set text(font: (raw-font, hei), size: 10pt)
  show link: underline
  show link: set text(blue)

  // 使用 i-figure 包添加章节编号
  //show heading: i-figured.reset-counters
  show figure: i-figured.show-figure
  show math.equation: i-figured.show-equation

  body
}


// 方程标注函数
#let pinit-highlight-equation-from(height: 2em, pos: bottom, fill: rgb(0, 180, 255), highlight-pins, point-pin, body) = {
  pinit-highlight(..highlight-pins, dy: -0.6em, fill: rgb(..fill.components().slice(0, -1), 40))
  pinit-point-from(
    fill: fill, pin-dx: -0.6em, pin-dy: if pos == bottom { 0.8em } else { -0.6em }, body-dx: 0pt, body-dy: if pos == bottom { -1.7em } else { -1.6em }, offset-dx: -0.6em, offset-dy: if pos == bottom { 0.8em + height } else { -0.6em - height },
    point-pin,
    rect(
      inset: 0.5em,
      stroke: (bottom: 0.12em + fill),
      {
        set text(fill: fill)
        body
      }
    )
  )
}
