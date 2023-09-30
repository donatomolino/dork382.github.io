function Notify(text, time) {

    // If a notification banner is already on the page, remove it
    var bannerEsistente = document.getElementById("notification-banner");
    if (bannerEsistente) {
      bannerEsistente.remove();
    }
    // Ceate the banner element 
    var banner = document.createElement("div");
    banner.id = "notification-banner";
    banner.classList.add("notification-banner");
    
    // Style the banner
    banner.style.position = "absolute";
    banner.style.top = "-100px";
    banner.style.left = "50%";
    banner.style.transform = "translateX(-50%)";
    banner.style.padding = "10px";
    banner.style.backgroundColor = "rgba(0, 0, 0, 0.7)";
    banner.style.color = "#fff";
    banner.style.borderRadius = "5px";
    banner.style.boxShadow = "0 0 10px rgba(0, 0, 0, 0.3)";
    banner.textContent = text;
  
    // Add the banner to the page
    document.body.appendChild(banner);
  
   // Animate the banner
  banner.style.transition = "top 0.5s ease-in-out";
  setTimeout(function() {
    banner.style.top = "10px"; 
  }, 100);

  // Hide the banner after a few seconds
  setTimeout(function() {
    banner.style.top = "-100px"; 
    setTimeout(function() {
      banner.remove(); 
    }, 500); 
  }, time * 1000);
}
  