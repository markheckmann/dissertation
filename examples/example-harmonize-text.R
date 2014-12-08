srt <- round(seq(0, 360, length=10))   # some angles

# plot text with standard rotation
plot(srt, rep(0, 10), type="n")
for (i in seq_along(srt))             
  text(srt[i], .5, label=srt[i], srt=srt[i])

# all labels in readable direction
hsrt <- harmonize_text_srt(srt)      
for (i in seq_along(srt))
  text(srt[i], 0, label=srt[i], srt=hsrt[i])

# all labels in readable direction but orthogonal to the ones before
hsrt <- harmonize_text_srt(srt, ortho = TRUE)  # harmonize rotation
for (i in seq_along(srt))
  text(srt[i],-.5, label=srt[i], srt=hsrt[i])

