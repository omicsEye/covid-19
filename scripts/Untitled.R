dt1 <- as.data.table(copy(state.x77))
dt1$state <- tolower(rownames(state.x77))
dt1 <- dt1[,.(state, Population)]
# only need state name and variable to plot in the input file:
str(dt1) 
us_map <- usmap::us_map() # used to add map scale

usmap::plot_usmap(data = dt1, values = "Population", labels = T)+
  labs(fill = 'State Population (1975)') + 
  scale_fill_gradientn(colours=rev(heat.colors(10)),na.value="grey90",
                       guide = guide_colourbar(barwidth = 25, barheight = 0.4,
                                               #put legend title on top of legend
                                               title.position = "top")) +
  # map scale
  ggsn::scalebar(data = us_map, dist = 500, dist_unit = "km",
                 border.size = 0.4, st.size = 4,
                 box.fill = c('black','white'),
                 transform = FALSE, model = "WGS84") + 
  # put legend at the bottom, adjust legend title and text font sizes
  theme(legend.position = "bottom",
        legend.title=element_text(size=12), 
        legend.text=element_text(size=10))

dt2 <- as.data.table(copy(state.x77))
dt2$state <- tolower(rownames(state.x77))
dt2 <- dt2[,.(state, Population)]
setkey(dt2, state)

states <- setDT(ggplot2::map_data("state"))
setkey(states, region)
# join data to map: left join states to dt2
dt2 <- dt2[states]
# data look like this: 
rmarkdown::paged_table(dt2[1:500,])

# create states location and abbreviations for label
# incl `Population` (the value to plot) in the label dataset, if want to fill with color. 
state_label_dt <- unique(dt2[, .(Population, x = mean(range(long)), y = mean(range(lat))), by = state])
snames <- data.table(state = tolower(state.name), abb = state.abb) # these are dataset within R
setkey(state_label_dt, state)
setkey(snames, state)
state_label_dt <- snames[state_label_dt]


# All labels for states to the right of lon = -77 will be on the right of lon = -50.
x_boundary = -77
x_limits <- c(-50, NA) # optional, for label repelling

ggplot(data = dt2, aes(x=long, y=lat, group=group))+
  geom_polygon(aes(fill=Population))+
  geom_path()+ 
  scale_fill_gradientn(colours = rev(heat.colors(10)),na.value = "grey90",
                       guide = guide_colourbar(barwidth = 25, barheight = 0.4,
                                               #put legend title on top of legend
                                               title.position = "top")) +
  # if need to repel labels... could further finetune
  geom_label_repel(data = state_label_dt[x>=x_boundary,],
                   aes(x = x,y = y, label = abb, fill = Population),
                   arrow = arrow(length = unit(0.02, "npc"), ends = "first"),
                   force = 5, hjust = 1, size = 3,
                   xlim  = x_limits, inherit.aes = F
  ) +
  # the normal labels: 
  geom_text(data=state_label_dt[x<x_boundary,], aes(x=x,y=y, label=abb), 
            size=3, inherit.aes=F) +
  coord_map() + 
  theme_classic() + 
  labs(fill = "Population", x = "Longitude", y = "Latitude") + 
  
  # map scale
  ggsn::scalebar(data = dt2, dist = 500, dist_unit = "km",
                 border.size = 0.4, st.size = 4,
                 box.fill = c('black','white'),
                 transform = TRUE, model = "WGS84") + 
  # put legend at the bottom, adjust legend title and text font sizes
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),  # font size of the legend 
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),  # remove axis, title, ticks
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line=element_blank())