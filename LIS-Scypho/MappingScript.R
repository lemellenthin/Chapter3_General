#Mapping

#Tutorials

#link1
#  https://www.molecularecologist.com/2016/07/01/making-maps-in-r-volume-2-ggplots/
# notes

# maptype
  #"terrain-background","watercolor","toner-lite")
base = get_map(location=c(-74,40.5,-71,41.4), zoom=10, maptype="terrain-background", source = c("google", "osm", "stamen"))
map1 = ggmap(base)
map1

# points on map individually
map1 + geom_point(data=cy2023, aes(x=longitude, y=latitude), 
                  color="white", cex=2.5) + # plot the points
  geom_point(data = visitnosight, aes(x=longitude, y=latitude), color="red", cex=2.5) +
  geom_point(data = visitsight, aes(x=longitude, y=latitude), color="green", cex=2.5) +
  geom_point(data = visitcollect, aes(x=longitude, y=latitude), color="blue", cex=2.5) +
  geom_point(data = visitDNA, aes(x=longitude, y=latitude), color="orange", cex=2.5) +
  geom_point(data = visitsequence, aes(x=longitude, y=latitude), color="pink", cex=2.5) +
 labs(x="Latitude", y="Longitude", title="Collection sites") + # label the axes
  theme_bw() + theme(legend.position="bottom", 
                     axis.text = element_text(size = rel(0.75)), 
                     legend.key = element_rect(colour = "white"), 
                     axis.text.x = element_text(angle=45, vjust=0.5)) # tweak the plot's appearance and legend position
#saving image
#pdf(file = "./Figures/1.pdf",   # The directory you want to save the file in
 #   width = 10, # The width of the plot in inches
  #  height = 6) # The height of the plot in inches

# points on map by group
map1 +
  geom_point(data=cy2023,aes(x=longitude,y=latitude,colour = Collections, 
                             shape = Collections, size=Collections)) + 
  scale_color_manual(values=c("red", "blue", "purple", "orange","green"),
                     labels=c("Visited Not Sighted (n=12)","Sighted (n=13)",
                              "Collected (n=6)",
                              "Collected, DNA (n=5)",
                              "Collected, DNA & Sequenced (n=1)")) + 
  scale_shape_manual(values=c(7,10,18,15,19),
                     labels=c("Visited Not Sighted (n=12)","Sighted (n=13)",
                             "Collected (n=6)",
                             "Collected, DNA (n=5)",
                             "Collected, DNA & Sequenced (n=1)")) + 
  scale_size_manual(values=c(5,5,3,4,9),
                    labels=c("Visited Not Sighted (n=12)","Sighted (n=13)",
                            "Collected (n=6)",
                            "Collected, DNA (n=5)",
                            "Collected, DNA & Sequenced (n=1)")) +
  labs(x="Latitude", y="Longitude", title="Collection sites 2023 \n(n=25)") +
  theme(plot.title=element_text(hjust=0.5, size=rel(2))) +
  theme(legend.title.align=0.5,
        panel.border = element_rect(colour = "black", fill=NA),
        legend.box.background = element_rect(color = "black",fill = "white"),
        legend.position=c(0.88, 0.292), 
        axis.text = element_text(size = rel(0.9)), 
        legend.key = element_rect(colour = "white"), 
        axis.text.x = element_text(angle=45, vjust=0.5),
        legend.text=element_text(size=8.45),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(1.4, 'cm'),
        legend.title=element_text(size=11))

#saving image
#dev.off()


#tried to add date labels
  geom_text_repel(data = cy2023, 
                  box.padding = 0.5
                  ,aes(x = longitude, y = latitude, 
                                     label = date, angle=90,vjust=0))
            #,size = 3, vjust = 0, hjust = -0.5, angle=20)

cy2023$date <- as.factor(cy2023$date)

#link2
  

# get all occurrence
spnames <- c('Cyanea capillata')
df_mult <- occ(query = spnames, from = c('gbif', 'ecoengine', 'bison', 'vertnet', 'inat'),
               limit = 100, geometry = geom, has_coords = T) # change limit to 10,000

df_nam <- fixnames(df_mult, how = "query")
df_comb <- occ2df(df_nam)

records_geoclean <- df_comb %>% coord_incomplete(drop = T) %>% coord_unlikely(drop = T)
# Want to do coord_impossible too, but some unkown error is happening
# records_geoclean <- df_comb %>% coord_impossible(drop = T) %>% coord_incomplete(drop = T) %>% coord_unlikely(drop = T)

records_dateclean <- records_geoclean
records_dateclean <- date_standardize(records_dateclean, "%Y")
records_dateclean <- date_missing(records_dateclean, format = "%Y", drop = T)
records_dateclean$date <- as.numeric(records_dateclean$date)

# Set which years for which you want data
records <- subset(records_dateclean, records_dateclean$date >= 1960 & records_dateclean$date <= 2017) # for modeling, make years go from 1960-1990

nrow(records)

records$longitude <- as.numeric(records$longitude)
records$latitude <- as.numeric(records$latitude)

recordsSpatial <- SpatialPointsDataFrame(coords = cbind(records$longitude, records$latitude), data = records,
                                         proj4string = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) # Turns our coordinates into an appropriate data frame

saveRDS(recordsSpatial, "/Test/SpatialRecords")
##  